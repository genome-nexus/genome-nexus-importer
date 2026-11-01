#!/usr/bin/env python3
"""Convert HGVS genomic notation to Ensembl default format for offline VEP.

VEP's cache/offline mode cannot read HGVS as input -- the Ensembl documentation
lists "Using HGVS notation as input (--format hgvs)" among the things the cache
does not support, and a FASTA file only enables HGVS as *output*.  To annotate a
set of HGVSg keys offline they must first be rewritten as Ensembl default format:

    chrom  start  end  ref/alt  strand  identifier

Two properties of that format carry the work:

  * ``ref`` is written as the placeholder ``N`` and resolved by VEP's
    ``--lookup_ref`` against the FASTA, so no reference sequence access is
    needed here.
  * ``identifier`` (the 6th field) is echoed back by VEP as the ``id`` of the
    annotation, so passing the *original, unmodified* HGVSg string makes the
    resulting document's ``_id`` match the source key exactly.

Deliberate differences from ``VEPController.hgvsgToEnsembl()`` in
genome-nexus-vep, which this replaces:

  * A deletion that states its deleted sequence (``g.100_102delGGG``) has that
    sequence **dropped**; the coordinate range is authoritative.  The Java
    version assigns it to the *alt* allele, producing ``N/GGG`` which
    ``--lookup_ref`` resolves to ``GGG/GGG`` -- a null variant rather than a
    deletion.  Verified against gn-vep-111-grch38: ``del``, ``delAGT`` and even
    a deliberately wrong ``delTTT`` all return byte-identical annotations, so
    dropping the sequence is lossless.
  * The compound ``g.100_102delGGGinsAT`` form is handled as a delins.  The Java
    regex rejects it outright (its trailing group is ``[A-Z]*``, which cannot
    match the lowercase ``ins``).
  * Inputs that genuinely need reference sequence (``dup``, ``inv``) or are
    malformed are written to a deferral file rather than silently dropped, so a
    later database-mode pass can pick them up.  The Java version throws, and the
    HTTP layer discards the error whenever at least one variant in the batch
    succeeded.

Deletions, insertions and substitutions are already minimal as emitted -- verified
byte-identical to database-mode output across 400 real production keys. ``delins``
is not: writing it over the full HGVS range leaves a redundant shared prefix or
suffix (``1:g.100011882_100011883delinsTT`` becomes ``CT/TT`` where database mode
reports ``C/T``). The annotation is equivalent either way -- hgvsc, hgvsp and the
consequence terms all match -- but ``allele_string`` and ``variant_class`` differ
from production. Supplying ``--fasta`` resolves the reference and trims to the
minimal form; without it, ``delins`` is emitted un-trimmed.

Output lines are written in input order.  Nothing is deduplicated: MongoDB
``_id`` values are unique by definition, so a key list projected from ``_id``
cannot contain duplicates.
"""

import argparse
import re
import sys
from collections import Counter

# chrom is deliberately permissive: production keys include 1-22, X, Y, M, MT
# and assembly contigs such as GL000205.2 / KI270438.1.  Anything VEP does not
# recognise fails at annotation time and is picked up as a straggler, which is
# preferable to guessing here.
_HGVSG = re.compile(r"^(?P<chrom>.+?):g\.(?P<start>\d+)(?:_(?P<end>\d+))?(?P<change>.*)$")

# Order matters only in that the ins-bearing delins pattern is tried before the
# plain deletion pattern.  The sequence classes are uppercase-only, so they
# cannot swallow the lowercase "ins"/"del"/"dup"/"inv" keywords.
_SUB    = re.compile(r"^(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$")
_DELINS = re.compile(r"^del(?P<delseq>[ACGTN]*)ins(?P<alt>[ACGTN]+)$")
_DEL    = re.compile(r"^del(?P<delseq>[ACGTN]*)$")
_INS    = re.compile(r"^ins(?P<alt>[ACGTN]+)$")
_DEFER  = re.compile(r"^(dup|inv)")

STRAND = "+"
REF_PLACEHOLDER = "N"   # resolved by --lookup_ref


class Faidx(object):
    """Minimal random-access reader for an uncompressed FASTA with a .fai index.

    Deliberately dependency-free: the VEP image carries no python3 at all, let
    alone pysam or pyfaidx, so the backfill image installs python3 and this
    reads the plain FASTA directly. A bgzip-compressed FASTA is not supported --
    decompress it once during cache staging.
    """

    def __init__(self, path):
        self.path = path
        self.index = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fh = open(path, "rb")

    def fetch(self, chrom, start, end):
        """1-based inclusive. Returns uppercase bases, or None if out of range."""
        rec = self.index.get(chrom)
        if rec is None:
            return None
        length, offset, linebases, linewidth = rec
        if start < 1 or end > length or end < start:
            return None
        # .fai gives the byte offset of the sequence; newlines must be skipped.
        def byte_at(pos):
            return offset + (pos - 1) // linebases * linewidth + (pos - 1) % linebases
        self._fh.seek(byte_at(start))
        raw = self._fh.read(byte_at(end) - byte_at(start) + 1)
        return raw.decode("ascii", "replace").replace("\n", "").replace("\r", "").upper()


def minimise(start, end, ref, alt):
    """Trim a shared suffix then prefix, the way VEP's HGVS parser does.

    Reproduces database-mode output exactly on the known divergent cases:
      100011882-100011883 CT/TT   -> 100011882-100011882 C/T
      100022945-100022946 GA/A    -> 100022945-100022945 G/-
      100015300-100015303 AGGT/TTT-> 100015300-100015302 AGG/TT
    """
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
        end -= 1
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        start += 1
    if not ref:
        # Became an insertion: Ensembl format places those with start = end + 1.
        return end + 1, end, "-", (alt or "-")
    if not alt:
        return start, end, ref, "-"
    return start, end, ref, alt


class Deferred(Exception):
    """Raised for input that must go to the database-mode pass instead."""


def convert(hgvsg, fasta=None):
    """Return an Ensembl default-format line for *hgvsg*.

    *fasta* is an optional :class:`Faidx`; when given, ``delins`` is trimmed to
    its minimal representation so ``allele_string`` matches database mode.

    Raises :class:`Deferred` for input that needs reference sequence access or
    that cannot be parsed.  The caller is expected to route those to a deferral
    file rather than drop them.
    """
    m = _HGVSG.match(hgvsg)
    if not m:
        raise Deferred("unparseable")

    chrom = m.group("chrom")
    start = int(m.group("start"))
    end = int(m.group("end")) if m.group("end") else start
    change = m.group("change")

    if end < start:
        raise Deferred("end before start")

    if _DEFER.match(change):
        # dup needs the duplicated sequence, inv needs it reverse-complemented.
        raise Deferred(change[:3])

    sub = _SUB.match(change)
    if sub:
        return _line(chrom, start, end, REF_PLACEHOLDER, sub.group("alt"), hgvsg)

    delins = _DELINS.match(change)
    if delins:
        # The deleted sequence is redundant with the range; the inserted one is
        # the variant and must be kept. Trimming needs the real reference, so it
        # only happens when a FASTA is supplied.
        alt = delins.group("alt")
        if fasta is not None:
            ref = fasta.fetch(chrom, start, end)
            if ref is not None and len(ref) == end - start + 1:
                m_start, m_end, m_ref, m_alt = minimise(start, end, ref, alt)
                return _line(chrom, m_start, m_end, m_ref, m_alt, hgvsg)
        return _line(chrom, start, end, REF_PLACEHOLDER, alt, hgvsg)

    dele = _DEL.match(change)
    if dele:
        # Stated deleted sequence intentionally ignored -- see module docstring.
        return _line(chrom, start, end, REF_PLACEHOLDER, "-", hgvsg)

    ins = _INS.match(change)
    if ins:
        # Ensembl format places an insertion between two bases by giving a start
        # one greater than the end, so the HGVS coordinates are swapped.  HGVS
        # requires the two flanking positions to be adjacent; anything else is
        # malformed and cannot be placed unambiguously.
        if end != start + 1:
            raise Deferred("insertion not between adjacent bases")
        return _line(chrom, end, start, "-", ins.group("alt"), hgvsg)

    raise Deferred("unrecognised change: %s" % (change[:20] or "empty"))


def _line(chrom, start, end, ref, alt, identifier):
    return "%s %d %d %s/%s %s %s" % (chrom, start, end, ref, alt, STRAND, identifier)


def run(infile, outfile, deferfile, fasta=None):
    stats = Counter()
    reasons = Counter()
    for raw in infile:
        hgvsg = raw.strip()
        if not hgvsg:
            continue
        stats["read"] += 1
        try:
            outfile.write(convert(hgvsg, fasta) + "\n")
            stats["converted"] += 1
        except Deferred as exc:
            deferfile.write(hgvsg + "\n")
            stats["deferred"] += 1
            reasons[str(exc)] += 1
    return stats, reasons


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="infile", required=True,
                    help="file of HGVSg keys, one per line ('-' for stdin)")
    ap.add_argument("--out", dest="outfile", required=True,
                    help="Ensembl default-format output for offline VEP ('-' for stdout)")
    ap.add_argument("--defer", dest="deferfile", required=True,
                    help="HGVSg keys needing the database-mode pass")
    ap.add_argument("--fasta", dest="fasta", default=None,
                    help="uncompressed reference FASTA with a .fai index; enables "
                         "minimal-representation trimming of delins")
    args = ap.parse_args(argv)

    def _open(path, mode):
        if path == "-":
            return sys.stdin if "r" in mode else sys.stdout
        return open(path, mode)

    fasta = Faidx(args.fasta) if args.fasta else None
    fin = _open(args.infile, "r")
    fout = _open(args.outfile, "w")
    fdefer = _open(args.deferfile, "w")
    try:
        stats, reasons = run(fin, fout, fdefer, fasta)
    finally:
        for f in (fin, fout, fdefer):
            if f not in (sys.stdin, sys.stdout):
                f.close()

    # Accounting goes to stderr so stdout stays pipeable into VEP.
    sys.stderr.write("read %d  converted %d  deferred %d%s\n"
                     % (stats["read"], stats["converted"], stats["deferred"],
                        "" if fasta else "  (no --fasta: delins left un-trimmed)"))
    for reason, n in reasons.most_common():
        sys.stderr.write("  deferred[%s] %d\n" % (reason, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
