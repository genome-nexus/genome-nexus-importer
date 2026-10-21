"""Shared sequence-alignment utilities for translating a cancer-hotspot residue/position
from one transcript's numbering to another's.

Used by both:
  - overwrite_transcript_and_residue_to_mskcc_isoform.py (GRCh37 native transcript -> MSKCC's
    preferred GRCh37 transcript)
  - liftover_hotspots_to_grch38.py (GRCh37 -> GRCh38)

Both steps are really the same problem: given a hotspot defined on transcript A, find the
correct residue/position on transcript B (which may be a different isoform of the same gene,
or the same gene on a different genome build). Never assume a raw position carries over
unchanged - always verify via alignment, since two transcripts can differ by an internal
indel that shifts positions differently on each side of it.
"""

import gzip

from Bio import Align
from Bio.Align import substitution_matrices


def load_fasta_sequences(path):
    """Returns (seqs, versions): both dicts keyed by unversioned transcript stable ID."""
    seqs, versions = {}, {}
    with gzip.open(path, 'rt') as f:
        tid, version, chunks = None, None, []
        for line in f:
            if line.startswith('>'):
                if tid is not None:
                    seqs[tid] = ''.join(chunks)
                    versions[tid] = version
                tid, version = None, None
                for field in line.split():
                    if field.startswith('transcript:'):
                        versioned = field.split(':', 1)[1]
                        tid, _, version = versioned.partition('.')
                chunks = []
            else:
                chunks.append(line.strip())
        if tid is not None:
            seqs[tid] = ''.join(chunks)
            versions[tid] = version
    return seqs, versions


_aligner = Align.PairwiseAligner()
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.open_gap_score = -11
_aligner.extend_gap_score = -1
_aligner.mode = "global"


def build_position_map(seq_a, seq_b):
    """1-indexed position in seq_a -> 1-indexed position in seq_b, or None if gapped."""
    alignment = _aligner.align(seq_a, seq_b)[0]
    aligned_a, aligned_b = alignment[0], alignment[1]
    pos_map = {}
    a_idx, b_idx = 0, 0
    for ca, cb in zip(aligned_a, aligned_b):
        if ca != '-':
            a_idx += 1
        if cb != '-':
            b_idx += 1
        if ca != '-':
            pos_map[a_idx] = b_idx if cb != '-' else None
    return pos_map


def map_position(pos_map, pos, seq_b):
    new_pos = pos_map.get(pos)
    if new_pos is None:
        return None, None
    new_aa = seq_b[new_pos - 1] if 1 <= new_pos <= len(seq_b) else None
    return new_pos, new_aa


def translate_residue(residue, rtype, pos_map, seq_from, seq_to):
    """Translates one hotspot's residue/position from seq_from's numbering to seq_to's.

    Returns (fields_or_None, status, detail). fields_or_None is {'residue': ..,
    'amino_acid_position': ..} on success, else None. status is one of:
    'unchanged' | 'shifted' | 'deleted' | 'grch37_ref_mismatch' | 'aa_changed' | 'unknown_type'.
    (the 'grch37_ref_mismatch' name is historical - it just means "seq_from doesn't have the
    stated reference amino acid at the stated position", regardless of which build seq_from is from)
    """
    if rtype in ("single residue", "3d"):
        letter, pos = residue[0], int(residue[1:])
        old_aa = seq_from[pos - 1] if 1 <= pos <= len(seq_from) else None
        if old_aa is None or old_aa != letter:
            return None, "grch37_ref_mismatch", f"seq_from has {old_aa} at {pos}, hotspot says {letter}"
        new_pos, new_aa = map_position(pos_map, pos, seq_to)
        if new_pos is None:
            return None, "deleted", f"position {pos} not present in target transcript"
        if new_aa != letter:
            return None, "aa_changed", f"position {pos}->{new_pos}, but aa is {new_aa} not {letter}"
        fields = {"residue": f"{letter}{new_pos}", "amino_acid_position": str(new_pos)}
        return fields, ("unchanged" if new_pos == pos else "shifted"), f"{pos}->{new_pos}"

    if rtype == "splice site":
        letter, pos = residue[0], int(residue[1:])
        new_pos, _new_aa = map_position(pos_map, pos, seq_to)
        if new_pos is None:
            return None, "deleted", f"position {pos} not present in target transcript"
        fields = {"residue": f"{letter}{new_pos}", "amino_acid_position": str(new_pos)}
        return fields, ("unchanged" if new_pos == pos else "shifted"), f"{pos}->{new_pos}"

    if rtype == "in-frame indel":
        if "-" in residue:
            start, end = (int(x) for x in residue.split("-"))
        else:
            start = end = int(residue)
        new_start, _ = map_position(pos_map, start, seq_to)
        new_end, _ = map_position(pos_map, end, seq_to)
        if new_start is None or new_end is None:
            return None, "deleted", f"range {residue} partially or fully absent in target transcript"
        new_residue = str(new_start) if new_start == new_end else f"{new_start}-{new_end}"
        fields = {"residue": new_residue, "amino_acid_position": new_residue}
        return fields, ("unchanged" if (new_start, new_end) == (start, end) else "shifted"), f"{residue}->{new_residue}"

    return None, "unknown_type", f"unhandled type {rtype}"
