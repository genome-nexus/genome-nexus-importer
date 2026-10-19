#!/usr/bin/env python3
"""Lift over a combined GRCh37 2D/3D cancer-hotspots file to GRCh38.

Replaces the previous approach (update_hotspots_to_grch38.py), which fetched whole
translated protein sequences from the Ensembl REST API and dropped an entire gene's
hotspots if the two full-length sequences weren't byte-identical anywhere.

Instead, this does a real per-residue liftover: for every gene, align the GRCh37 and
GRCh38 protein sequences (global alignment, BLOSUM62) once per transcript pair, then walk
that alignment to find exactly where each hotspot residue lands in GRCh38 - whether
that's the same position, a shifted position (including a different shift on each side
of an internal exon difference), or no position at all (the residue's exon isn't present
in the GRCh38 transcript). It also detects the case where the raw numeric position looks
unchanged but the two transcripts actually diverge elsewhere, which the alignment can
tell apart from a position that's merely sitting in a region of true 1:1 identity.

A single fixed transcript per gene isn't always right on either side, so both sides try
multiple candidates and pick whichever actually reproduces the hotspot's stated reference
amino acid:
  - GRCh37 side: the transcript currently on the row (after the MSKCC GRCh37 isoform
    override step), then that row's original pre-override transcript (cancerhotspots.org's
    own choice) as a fallback - covers cases like a 3D-sourced hotspot whose native
    transcript differs from the gene's general 2D override.
  - GRCh38 side: every MSKCC override entry for the gene (a gene can have more than one
    row, e.g. CDKN2A's two reading frames), then MANE Select as a final fallback.

The actual position-translation logic (build an alignment, walk it, check the hotspot's
reference amino acid survives the trip) lives in hotspot_transcript_utils.py and is shared
with overwrite_transcript_and_residue_to_mskcc_isoform.py, which has the identical problem
one step earlier (GRCh37 native transcript -> MSKCC's preferred GRCh37 transcript).
"""

import argparse
import csv
import sys
from collections import Counter

from hotspot_transcript_utils import load_fasta_sequences, build_position_map, translate_residue


def load_grch38_overrides(path):
    """gene_name -> ordered list of candidate (unversioned) transcript IDs.

    A gene may legitimately have more than one row (e.g. CDKN2A's two reading frames) -
    keep all of them, in file order, rather than collapsing to one.
    """
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["enst_id"]:
                out.setdefault(row["gene_name"], []).append(row["enst_id"].split(".")[0])
    return out


def load_mane_select(path):
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["MANE_status"] == "MANE Select":
                out[row["symbol"]] = row["Ensembl_nuc"].split(".")[0]
    return out


def load_hgnc_alias_map(path):
    current_symbols, prev_to_current = set(), {}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("symbol") or ""
            current_symbols.add(sym)
            for field in ("prev_symbol", "alias_symbol"):
                for p in (row.get(field) or "").split("|"):
                    p = p.strip()
                    if p:
                        prev_to_current.setdefault(p, sym)
    return current_symbols, prev_to_current


def load_original_grch37_transcripts(path):
    """(gene, residue) -> cancerhotspots.org's own pre-MSKCC-override transcript ID."""
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[(row["hugo_symbol"], row["residue"])] = row["transcript_id"]
    return out


def find_grch38_candidates(gene, overrides, mane_select, current_symbols, prev_to_current,
                            mane_select_casefold):
    lookup_gene = gene if gene in current_symbols else prev_to_current.get(gene, gene)
    candidates = []
    for tid in overrides.get(lookup_gene, []):
        candidates.append((tid, "mskcc_override"))
    if lookup_gene in mane_select:
        candidates.append((mane_select[lookup_gene], "mane_select"))
    elif lookup_gene.casefold() in mane_select_casefold:
        # a few genes (the "C#orfN" family) only differ by case between data sources
        # (e.g. hotspots' "C10ORF71" vs MANE's "C10orf71")
        candidates.append((mane_select_casefold[lookup_gene.casefold()], "mane_select_casefold"))
    return candidates


def liftover_row_trying_all_candidates(row, grch37_candidates, grch38_candidates,
                                        grch37_seqs, grch38_seqs, grch38_versions,
                                        pos_map_cache):
    """Tries every (grch37_candidate, grch38_candidate) combination in order and returns
    the first one that succeeds. If none succeed, returns the first attempt's failure
    (generally the most informative: current transcript on current preferred target).
    """
    residue, rtype = row["residue"], row["type"]
    first_failure = None
    for tid37 in grch37_candidates:
        seq37 = grch37_seqs.get(tid37)
        if seq37 is None:
            if first_failure is None:
                first_failure = (None, "no_sequence", f"no grch37 sequence for {tid37}")
            continue
        for tid38, _source in grch38_candidates:
            seq38 = grch38_seqs.get(tid38)
            if seq38 is None:
                if first_failure is None:
                    first_failure = (None, "no_sequence", f"no grch38 sequence for {tid38}")
                continue

            cache_key = (tid37, tid38)
            if cache_key not in pos_map_cache:
                pos_map_cache[cache_key] = build_position_map(seq37, seq38)
            pos_map = pos_map_cache[cache_key]

            fields, status, detail = translate_residue(residue, rtype, pos_map, seq37, seq38)
            if fields is not None:
                lifted = dict(row)
                lifted.update(fields)
                lifted["transcript_id"] = tid38
                lifted["transcript_id_version"] = grch38_versions.get(tid38)
                return lifted, status, detail
            if first_failure is None:
                first_failure = (None, status, detail)
    return first_failure if first_failure else (None, "no_grch38_transcript", "no candidates available")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--grch37-hotspots", required=True, help="grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt")
    parser.add_argument("--grch37-hotspots-original", required=True,
                        help="common_input/hotspots/hotspots_v2_v3_and_3d_grch37.txt (pre-MSKCC-override, cancerhotspots.org's own transcript choices)")
    parser.add_argument("--grch37-pep-fasta", required=True, help="release-pinned Homo_sapiens.GRCh37.pep.all.fa.gz")
    parser.add_argument("--grch38-pep-fasta", required=True, help="release-pinned Homo_sapiens.GRCh38.pep.all.fa.gz")
    parser.add_argument("--grch38-overrides", required=True, help="isoform_overrides_at_mskcc_grch38.txt")
    parser.add_argument("--mane-select", required=True, help="MANE.GRCh38.vX.Y.summary.txt")
    parser.add_argument("--hgnc", required=True, help="hgnc_complete_set file, for symbol normalization fallback")
    parser.add_argument("--out", required=True, help="output path for the lifted-over grch38 hotspots file")
    parser.add_argument("--dropped-out", required=True, help="output path for rows that could not be lifted over, with reasons")
    args = parser.parse_args()

    print("Loading sequences and lookup tables...", file=sys.stderr)
    grch37_seqs, _grch37_versions = load_fasta_sequences(args.grch37_pep_fasta)
    grch38_seqs, grch38_versions = load_fasta_sequences(args.grch38_pep_fasta)
    overrides_38 = load_grch38_overrides(args.grch38_overrides)
    mane_select = load_mane_select(args.mane_select)
    mane_select_casefold = {k.casefold(): v for k, v in mane_select.items()}
    current_symbols, prev_to_current = load_hgnc_alias_map(args.hgnc)
    original_transcripts = load_original_grch37_transcripts(args.grch37_hotspots_original)

    with open(args.grch37_hotspots, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)

    pos_map_cache = {}
    gene_grch38_candidates = {}

    kept_rows = []
    dropped_rows = []

    for row in rows:
        gene = row["hugo_symbol"]
        tid37_current = row["transcript_id"]
        tid37_original = original_transcripts.get((gene, row["residue"]))

        grch37_candidates = [tid37_current]
        if tid37_original and tid37_original != tid37_current:
            grch37_candidates.append(tid37_original)

        if gene not in gene_grch38_candidates:
            gene_grch38_candidates[gene] = find_grch38_candidates(
                gene, overrides_38, mane_select, current_symbols, prev_to_current, mane_select_casefold)
        grch38_candidates = gene_grch38_candidates[gene]

        if not grch38_candidates:
            dropped_rows.append({**row, "liftover_status": "no_grch38_transcript", "liftover_detail": "no mskcc override or MANE Select for this gene"})
            continue

        # A handful of cancerhotspots.org's own "GRCh37" transcript IDs don't exist in the
        # GRCh37 database at all - they're already GRCh38 transcript IDs (a pre-existing
        # quirk in the upstream data, not something introduced here). If so, there's nothing
        # to lift over: treat it as already-GRCh38 and keep it as-is.
        if all(grch37_seqs.get(t) is None for t in grch37_candidates) and tid37_current in grch38_seqs:
            out = dict(row)
            out["transcript_id_version"] = grch38_versions.get(tid37_current)
            out["liftover_status"] = "already_grch38"
            out["liftover_detail"] = f"transcript_id {tid37_current} is already a valid GRCh38 transcript"
            kept_rows.append(out)
            continue

        lifted, status, detail = liftover_row_trying_all_candidates(
            row, grch37_candidates, grch38_candidates, grch37_seqs, grch38_seqs, grch38_versions, pos_map_cache)

        if lifted is None:
            dropped_rows.append({**row, "liftover_status": status, "liftover_detail": detail})
        else:
            lifted["liftover_status"] = status
            lifted["liftover_detail"] = detail
            kept_rows.append(lifted)

    out_fields = fieldnames + ["liftover_status", "liftover_detail"]
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for r in kept_rows:
            writer.writerow(r)

    with open(args.dropped_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for r in dropped_rows:
            writer.writerow(r)

    print(f"Total rows: {len(rows)}", file=sys.stderr)
    print(f"  Kept: {len(kept_rows)}", file=sys.stderr)
    print(f"  Dropped: {len(dropped_rows)}", file=sys.stderr)
    print("  Dropped breakdown:", dict(Counter(r["liftover_status"] for r in dropped_rows)), file=sys.stderr)
    print("  Kept breakdown:", dict(Counter(r["liftover_status"] for r in kept_rows)), file=sys.stderr)
    print(f"Wrote {args.out} and {args.dropped_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
