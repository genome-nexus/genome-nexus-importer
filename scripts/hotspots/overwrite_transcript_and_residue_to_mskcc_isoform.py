#!/usr/bin/env python3
"""Swap in MSKCC's preferred GRCh37 transcript for each cancer hotspot gene, translating
the residue/position via sequence alignment rather than assuming the raw number carries
over unchanged.

A gene can have more than one candidate override transcript (e.g. CDKN2A's two reading
frames) - every candidate is tried, and whichever one actually reproduces the hotspot's
reference amino acid (at the same position, or a different one if the alignment finds a
shift) is used. If no candidate works, the row is left on its original transcript rather
than being silently corrupted - and reported separately so it gets reviewed, not missed.

The actual alignment/translation logic lives in hotspot_transcript_utils.py and is shared
with liftover_hotspots_to_grch38.py, which has the identical problem one step later
(GRCh37 -> GRCh38, instead of GRCh37 native transcript -> MSKCC's preferred one).
"""
import argparse
import sys
from collections import Counter

import pandas as pd

from hotspot_transcript_utils import load_fasta_sequences, build_position_map, translate_residue


def load_hgnc(hgnc_path: str):
    """
    Load HGNC file, return:
    - current_symbols: set of approved symbols
    - prev_to_current: map of old symbol -> current symbol
    """
    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str).fillna("")

    current_symbols = set(hgnc["symbol"].astype(str))

    prev_to_current = {}
    for _, row in hgnc.iterrows():
        current = row["symbol"]
        prev_field = row.get("prev_symbol", "")
        if not prev_field:
            continue
        for prev_sym in prev_field.split("|"):
            prev_sym = prev_sym.strip()
            if prev_sym:
                prev_to_current.setdefault(prev_sym, current)

    return current_symbols, prev_to_current


def normalize_symbol(symbol: str,
                     current_symbols: set,
                     prev_to_current: dict) -> str:
    """
    1. if symbol already current, return as-is
    2. else if symbol is a prev_symbol, map to current
    3. else keep as-is
    """
    if symbol in current_symbols:
        return symbol
    if symbol in prev_to_current:
        return prev_to_current[symbol]
    return symbol


def load_isoform_overrides(path: str) -> dict:
    """
    Read isoform_overrides_at_mskcc_grch37.txt and return:
      gene_name -> ordered list of candidate (unversioned) transcript IDs

    A gene may legitimately have more than one row (e.g. CDKN2A's two reading frames) -
    keep all of them, in file order, rather than collapsing to one.
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out = {}
    for _, row in df.iterrows():
        gene = row["gene_name"]
        enst_id = row["enst_id"]
        if enst_id:
            out.setdefault(gene, []).append(enst_id.split(".")[0])
    return out


def load_biomart_versions(biomart_path: str) -> dict:
    """
    Build:
      transcript_stable_id (no version) -> version number string
    from ensembl_biomart_geneids.txt
    """
    biomart = pd.read_csv(biomart_path, sep="\t", dtype=str).fillna("")
    versions = {}
    for _, row in biomart.iterrows():
        tid = row["Transcript stable ID"]
        ver_full = row["Versioned transcript ID"]
        if tid:
            if "." in ver_full:
                _, ver_num = ver_full.split(".", 1)
            else:
                ver_num = ""
            versions.setdefault(tid, ver_num)
    return versions


def load_residue_change_map(residue_path: str) -> dict:
    """
    Build nested dict:
      residue_change[transcript_id][original_residue] = revised_residue
    from mskcc_isoform_residue_changes.txt

    This is a manual escape hatch for edge cases the automatic alignment-based
    translation below can't resolve on its own - it should rarely be needed now, since
    alignment already finds the correct shift for the common case.
    """
    res_df = pd.read_csv(residue_path, sep="\t", dtype=str).fillna("")
    m = {}
    for _, row in res_df.iterrows():
        tid = row["transcript_id"]
        orig = row["original_residue"]
        revised = row["revised_residue"]
        if tid not in m:
            m[tid] = {}
        m[tid].setdefault(orig, revised)
    return m


def pick_override_transcript(row, candidates, seqs, versions, biomart_ver_map):
    """Tries each override candidate transcript in order and returns the first one whose
    sequence actually reproduces this hotspot's reference amino acid - either unchanged
    (identical sequence, or same position) or at a position found via alignment.

    Returns (new_tid, new_version, residue, amino_acid_position, status, detail).
    status is 'unchanged' | 'shifted' | <failure reason from translate_residue>.
    """
    orig_tid = str(row["transcript_id"])
    seq_from = seqs.get(orig_tid)

    for cand_tid in candidates:
        seq_to = seqs.get(cand_tid)
        if seq_to is None:
            continue
        cand_version = biomart_ver_map.get(cand_tid) or versions.get(cand_tid, "")

        if seq_from is not None and seq_to == seq_from:
            # identical protein - no translation needed, position carries over as-is
            return cand_tid, cand_version, row["residue"], row["amino_acid_position"], "unchanged", "identical sequence"

        if seq_from is None:
            continue

        pos_map = build_position_map(seq_from, seq_to)
        fields, status, detail = translate_residue(row["residue"], row["type"], pos_map, seq_from, seq_to)
        if fields is not None:
            return cand_tid, cand_version, fields["residue"], fields["amino_acid_position"], status, detail

    return None, None, None, None, "no_candidate_matched", f"none of {candidates} reproduced the reference amino acid"


def transform_hotspots(
    hotspots_df: pd.DataFrame,
    current_symbols: set,
    prev_to_current: dict,
    override_map: dict,
    biomart_ver_map: dict,
    residue_change_map: dict,
    seqs: dict,
    versions: dict,
) -> "tuple[pd.DataFrame, pd.DataFrame]":
    """
    For each row in hotspots_df:
      1. Normalize hugo_symbol using HGNC.
      2. Determine transcript_id + transcript_id_version by trying every MSK override
         candidate for the gene and keeping whichever one's sequence actually matches
         the hotspot (translating residue/position via alignment if needed). Falls back
         to the original transcript (unchanged) if no candidate matches.
      3. Apply any manual residue_change_map correction on top (escape hatch).

    Returns (transformed_df, review_df) - review_df lists every row where no override
    candidate could be verified, so it's reviewed rather than silently left wrong or
    silently dropped.
    """

    new_hugo_symbols = []
    new_transcript_ids = []
    new_versions = []
    new_residues = []
    new_positions = []
    review_rows = []

    status_counts = Counter()

    for _, row in hotspots_df.iterrows():
        original_symbol = str(row["hugo_symbol"])
        final_symbol = normalize_symbol(original_symbol, current_symbols, prev_to_current)

        orig_tid = str(row["transcript_id"])
        orig_residue = str(row["residue"])
        orig_position = str(row["amino_acid_position"])

        candidates = override_map.get(final_symbol, [])
        if candidates:
            final_tid, final_ver, final_residue, final_position, status, detail = pick_override_transcript(
                row, candidates, seqs, versions, biomart_ver_map)
            if final_tid is None:
                # no candidate verified - keep the original transcript rather than guess
                final_tid = orig_tid
                final_ver = biomart_ver_map.get(orig_tid, "")
                final_residue = orig_residue
                final_position = orig_position
                review_rows.append({**row, "review_reason": detail})
            status_counts[status] += 1
        else:
            final_tid = orig_tid
            final_ver = biomart_ver_map.get(orig_tid, "")
            final_residue = orig_residue
            final_position = orig_position
            status_counts["no_override"] += 1

        # manual escape hatch, applied on top - rarely needed now that translation is
        # automatic, but kept for any case that still needs a hand-curated correction
        if final_tid in residue_change_map and orig_residue in residue_change_map[final_tid]:
            revised = residue_change_map[final_tid][orig_residue]
            final_residue = revised
            digits = "".join(c for c in revised if c.isdigit() or c == "-")
            final_position = digits if digits else final_position

        new_hugo_symbols.append(final_symbol)
        new_transcript_ids.append(final_tid)
        new_versions.append(final_ver)
        new_residues.append(final_residue)
        new_positions.append(final_position)

    print("Transcript selection status counts:", dict(status_counts), file=sys.stderr)
    print(f"Rows needing review (no override candidate verified): {len(review_rows)}", file=sys.stderr)

    out = hotspots_df.copy()
    out["hugo_symbol"] = new_hugo_symbols
    out["transcript_id"] = new_transcript_ids
    out["transcript_id_version"] = new_versions
    out["residue"] = new_residues
    out["amino_acid_position"] = new_positions

    cols = list(out.columns)
    if "transcript_id_version" in cols and "transcript_id" in cols:
        cols.insert(cols.index("transcript_id") + 1, cols.pop(cols.index("transcript_id_version")))
        out = out[cols]

    review_df = pd.DataFrame(review_rows) if review_rows else pd.DataFrame(columns=list(hotspots_df.columns) + ["review_reason"])
    return out, review_df


def main():
    parser = argparse.ArgumentParser(
        description="Transform hotspots using HGNC symbol updates, "
                    "MSK isoform overrides (verified by sequence alignment), "
                    "Biomart transcript versions, and residue remaps."
    )
    parser.add_argument("--hgnc", required=True,
                        help="Path to hgnc_set.txt")
    parser.add_argument("--overrides", required=True,
                        help="Path to isoform_overrides_at_mskcc_grch37.txt")
    parser.add_argument("--residue_changes", required=True,
                        help="Path to mskcc_isoform_residue_changes.txt")
    parser.add_argument("--biomart", required=True,
                        help="Path to ensembl_biomart_gene_ids.txt")
    parser.add_argument("--pep-fasta", required=True,
                        help="release-pinned Homo_sapiens.GRCh37.pep.all.fa.gz, for verifying/translating override transcripts")
    parser.add_argument("--hotspots", required=True,
                        help="Path to hotspots.txt (input)")
    parser.add_argument("--out", required=True,
                        help="Path to write transformed hotspots TSV")
    parser.add_argument("--review-out", required=True,
                        help="Path to write rows where no override candidate could be verified (kept on original transcript)")

    args = parser.parse_args()

    # load inputs
    current_symbols, prev_to_current = load_hgnc(args.hgnc)
    override_map = load_isoform_overrides(args.overrides)
    biomart_ver_map = load_biomart_versions(args.biomart)
    residue_change_map = load_residue_change_map(args.residue_changes)
    seqs, versions = load_fasta_sequences(args.pep_fasta)
    hotspots_df = pd.read_csv(args.hotspots, sep="\t", dtype=str).fillna("")

    # transform
    result_df, review_df = transform_hotspots(
        hotspots_df,
        current_symbols,
        prev_to_current,
        override_map,
        biomart_ver_map,
        residue_change_map,
        seqs,
        versions,
    )

    # write output
    result_df.to_csv(args.out, sep="\t", index=False)
    review_df.to_csv(args.review_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
