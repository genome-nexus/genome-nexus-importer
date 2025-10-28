#!/usr/bin/env python3
import argparse
import pandas as pd


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
    Read mskcc_isoform_overrides.txt and return:
      gene_name -> enst_id (versioned)
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out = {}
    for _, row in df.iterrows():
        gene = row["gene_name"]
        enst_id = row["enst_id"]
        out.setdefault(gene, enst_id)
    return out


def split_version(enst_with_ver: str):
    """
    ENST00000257430.4 -> ("ENST00000257430", "4")
    ENST00000257430   -> ("ENST00000257430", "")
    """
    if "." in enst_with_ver:
        main, ver = enst_with_ver.split(".", 1)
        return main, ver
    else:
        return enst_with_ver, ""


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


def transform_hotspots(
    hotspots_df: pd.DataFrame,
    current_symbols: set,
    prev_to_current: dict,
    override_map: dict,
    biomart_ver_map: dict,
    residue_change_map: dict
) -> pd.DataFrame:
    """
    For each row in hotspots_df:
      1. Normalize hugo_symbol using HGNC.
      2. Determine transcript_id + transcirpt_id_version.
         - Prefer MSK override (gene_name -> enst_id).
         - Else keep original transcript_id and infer version via Ensembl.
      3. Possibly rewrite residue using residue_change_map.
      4. Add new column transcirpt_id_version (next to transcript_id).
    """

    new_hugo_symbols = []
    new_transcript_ids = []
    new_versions = []
    new_residues = []

    for _, row in hotspots_df.iterrows():
        # Step 1. normalize symbol
        original_symbol = str(row["hugo_symbol"])
        final_symbol = normalize_symbol(
            original_symbol,
            current_symbols,
            prev_to_current
        )

        # Step 2. transcript selection + version
        orig_tid = str(row["transcript_id"])

        if final_symbol in override_map:
            override_enst = override_map[final_symbol]  # e.g. ENSTxxx.y
            main_tid, ver_num = split_version(override_enst)
            final_tid = main_tid
            final_ver = ver_num
        else:
            final_tid = orig_tid
            final_ver = biomart_ver_map.get(orig_tid, "")

        # Step 3. residue remap
        orig_residue = str(row["residue"])
        final_residue = orig_residue
        if final_tid in residue_change_map:
            if orig_residue in residue_change_map[final_tid]:
                final_residue = residue_change_map[final_tid][orig_residue]

        new_hugo_symbols.append(final_symbol)
        new_transcript_ids.append(final_tid)
        new_versions.append(final_ver)
        new_residues.append(final_residue)

    out = hotspots_df.copy()
    out["hugo_symbol"] = new_hugo_symbols
    out["transcript_id"] = new_transcript_ids
    out["transcript_id_version"] = new_versions  # keep your spelling
    out["residue"] = new_residues

    cols = list(out.columns)
    if "transcript_id_version" in cols and "transcript_id" in cols:
        cols.insert(cols.index("transcript_id") + 1, cols.pop(cols.index("transcript_id_version")))
        out = out[cols]

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Transform hotspots using HGNC symbol updates, "
                    "MSK isoform overrides, Biomart transcript versions, "
                    "and residue remaps."
    )
    parser.add_argument("--hgnc", required=True,
                        help="Path to hgnc_set.txt")
    parser.add_argument("--overrides", required=True,
                        help="Path to mskcc_isoform_overrides.txt")
    parser.add_argument("--residue_changes", required=True,
                        help="Path to mskcc_isoform_residue_changes.txt")
    parser.add_argument("--biomart", required=True,
                        help="Path to ensembl_biomart_gene_ids.txt")
    parser.add_argument("--hotspots", required=True,
                        help="Path to hotspots.txt (input)")
    parser.add_argument("--out", required=True,
                        help="Path to write transformed hotspots TSV")

    args = parser.parse_args()

    # load inputs
    current_symbols, prev_to_current = load_hgnc(args.hgnc)
    override_map = load_isoform_overrides(args.overrides)
    biomart_ver_map = load_biomart_versions(args.biomart)
    residue_change_map = load_residue_change_map(args.residue_changes)
    hotspots_df = pd.read_csv(args.hotspots, sep="\t", dtype=str).fillna("")

    # transform
    result_df = transform_hotspots(
        hotspots_df,
        current_symbols,
        prev_to_current,
        override_map,
        biomart_ver_map,
        residue_change_map
    )

    # write output
    result_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
