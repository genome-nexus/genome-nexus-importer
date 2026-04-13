#!/usr/bin/env python3
"""Combine 2d and 3d hotspot files. Add mutation type counts per hotspot. Then filter out the hotspots not following the threshold."""

import argparse
import pandas as pd
import numpy as np
import sys
import re
import hgvs.parser
PARSER = hgvs.parser.Parser()

def count_variant_types(row):
    MISSENSE_INDEX = 0
    TRUNC_INDEX = 1
    # FS_INDEX = 2
    INFRAME_INDEX = 2
    SPLICE_INDEX = 3
    TOTAL_INDEX = 4
    INDEX_NAMES = "missense trunc inframe splice total".split()
    rv = [0, 0, 0, 0, 0]
    
    if row.type == "3d":
        # ignore 3d hotspot
        return pd.Series([np.nan]*len(INDEX_NAMES), index=INDEX_NAMES)

    for v in row.variant_amino_acid.split("|"):
        # add total
        rv[TOTAL_INDEX] += int(v.split(":")[1])
        
        if "*" in v:
            rv[TRUNC_INDEX] += int(v.split(":")[1])
        elif "sp" in v:
            rv[SPLICE_INDEX] += int(v.split(":")[1])
        else:
            # parse hgvs to get protein change length (can always use TP53, don't care about which protein)
            if row.type == "in-frame indel":
                sv = PARSER.parse_hgvs_variant("TP53:p.{}".format(v.split(":")[0]))
            elif row.type == "single residue":
                sv = PARSER.parse_hgvs_variant("TP53:p.{}{}".format(row.residue, v.split(":")[0]))
            else:
                raise(Exception("unknown Type: {}".format(row.type)))
        
            if sv.posedit.pos.start == sv.posedit.pos.end and sv.posedit.length_change() == 0:
                rv[MISSENSE_INDEX] += int(v.split(":")[1])
            else:
                rv[INFRAME_INDEX] += int(v.split(":")[1])
        
    return pd.Series(rv, index=INDEX_NAMES)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hotspots_2d", default="../data/hotspots/v2_multi_type_residue.txt", type=str, help="2D cancerhotspots data file (v2)")
    parser.add_argument("hotspots_3d", default="../data/hotspots/3d_hotspots.txt", type=str, help="3D cancerhotspots data file")
    parser.add_argument("--hotspots_2d_v3", default=None, type=str,
                        help="2D cancerhotspots v3 data file (v2+v3 combined). Rows not in hotspots_2d are tagged as v3.")
    parser.add_argument("--removed_hotspots", default=None, type=str, help='Output removed hotspots')
    parser.add_argument("--override_unassigned_transcript_id_2d_hotspots", default=None, required=True, type=str,
                        help='Override transcript_id field for 2d hotspots without assignment')
    parser.add_argument("--override_unassigned_transcript_id_3d_hotspots", default=None, required=True, type=str,
                        help='Override transcript_id field for 3d hotspots without assignment')
    args = parser.parse_args()


    hotspots_2d = pd.read_csv(args.hotspots_2d, sep="\t", dtype=str)
    hotspots_2d.columns = [c.lower().replace("-","_") for c in hotspots_2d.columns]
    hotspots_2d['type'] = hotspots_2d.indel_size.fillna(0).apply(lambda x: "in-frame indel" if int(x) > 0 else "single residue")
    hotspots_2d.loc[((hotspots_2d.type == "single residue") & hotspots_2d.residue.str.contains("X")), 'type'] = "splice site"
    hotspots_2d['version'] = 'v2'

    # If v3 data provided, identify v3-only rows and tag them
    if args.hotspots_2d_v3:
        hotspots_2d_v3 = pd.read_csv(args.hotspots_2d_v3, sep="\t", dtype=str)
        hotspots_2d_v3.columns = [c.lower().replace("-","_") for c in hotspots_2d_v3.columns]
        hotspots_2d_v3['type'] = hotspots_2d_v3.indel_size.fillna(0).apply(lambda x: "in-frame indel" if int(x) > 0 else "single residue")
        hotspots_2d_v3.loc[((hotspots_2d_v3.type == "single residue") & hotspots_2d_v3.residue.str.contains("X")), 'type'] = "splice site"

        # Identify v3-only rows: in v3 file but not in v2 file (by hugo_symbol + residue)
        v2_keys = set(zip(hotspots_2d.hugo_symbol, hotspots_2d.residue))
        is_v3_only = hotspots_2d_v3.apply(
            lambda row: (row.hugo_symbol, row.residue) not in v2_keys, axis=1
        )
        v3_only = hotspots_2d_v3[is_v3_only].copy()
        v3_only['version'] = 'v3'
        print(f"Found {len(v3_only)} v3-only hotspots", file=sys.stderr)

        # Inherit transcript_id from v2 hotspots for the same gene
        v2_gene_transcript = hotspots_2d.dropna(subset=['transcript_id']).drop_duplicates('hugo_symbol').set_index('hugo_symbol')['transcript_id']
        v3_missing_tid = v3_only.transcript_id.isna()
        v3_only.loc[v3_missing_tid, 'transcript_id'] = v3_only.loc[v3_missing_tid, 'hugo_symbol'].map(v2_gene_transcript)
        n_still_missing = v3_only.transcript_id.isna().sum()
        if n_still_missing > 0:
            print(f"WARNING: {n_still_missing} v3-only hotspots still missing transcript_id:", file=sys.stderr)
            print(v3_only[v3_only.transcript_id.isna()][['hugo_symbol','residue']].to_string(), file=sys.stderr)

        # Combine v2 + v3-only 2D hotspots
        hotspots_2d = pd.concat([hotspots_2d, v3_only])

    hotspots_3d = pd.read_csv(args.hotspots_3d, sep="\t", dtype=str)
    hotspots_3d.columns = [c.lower().replace("-","_") for c in hotspots_3d.columns]
    # add type column
    hotspots_3d['type'] = '3d'
    hotspots_3d['version'] = 'v2'

    hotspots = pd.concat([hotspots_2d, hotspots_3d])
    assert(len(hotspots) == len(hotspots_2d) + len(hotspots_3d))

    h_counts = hotspots.apply(count_variant_types, axis=1)
    for c in h_counts.columns:
        hotspots[c + "_count"] = h_counts[c]
    for c in h_counts.columns[:-1]:
        hotspots[c+"_fraction"] = h_counts[c] / h_counts["total"]

    # fill in transcript id for those without an assignment (uses mskcc transcript overrides)
    if args.override_unassigned_transcript_id_2d_hotspots:
        hotspots2_overrides = pd.read_csv(args.override_unassigned_transcript_id_2d_hotspots, sep="\t").set_index("hugo_symbol")
        hotspots2_without_transcript_id = ((hotspots.type != "3d") & pd.isnull(hotspots.transcript_id))
        hotspots.loc[hotspots2_without_transcript_id, 'transcript_id'] = hotspots[hotspots2_without_transcript_id]["hugo_symbol"].apply(lambda x: hotspots2_overrides.loc[x])

        # no unassigned hugo symbols should be left after assigning missing transcripts for hotspots v2
        # (uses uniprot overrides)
        hotspots2_without_transcript_id = ((hotspots.type != "3d") & pd.isnull(hotspots.transcript_id))
        assert(hotspots2_without_transcript_id.sum() == 0)

    # fill in transcript id for those without an assignment
    if args.override_unassigned_transcript_id_3d_hotspots:
        hotspots3d_overrides = pd.read_csv(args.override_unassigned_transcript_id_3d_hotspots, sep="\t").set_index("hugo_symbol")
        hotspots3d_without_transcript_id = ((hotspots.type == "3d") & pd.isnull(hotspots.transcript_id))
        hotspots.loc[hotspots3d_without_transcript_id, 'transcript_id'] = hotspots[hotspots3d_without_transcript_id]["hugo_symbol"].apply(lambda x: hotspots3d_overrides.loc[x])


        # no unassigned hugo symbols should be left after assigning missing transcripts for hotspots v2
        hotspots3d_without_transcript_id = ((hotspots.type != "3d") & pd.isnull(hotspots.transcript_id))
        assert(hotspots3d_without_transcript_id.sum() == 0)
    
    # Add transcript_id_version column (empty — populated downstream by isoform override scripts)
    if 'transcript_id_version' not in hotspots.columns:
        hotspots['transcript_id_version'] = np.nan

    # Ensure consistent column order matching the expected export schema
    output_cols = [
        'hugo_symbol', 'residue', 'reference_amino_acid', 'amino_acid_position',
        'cluster', 'pdb_chains', 'class', 'variant_amino_acid',
        'q_value', 'qvalue_pancan', 'qvaluect', 'p_value',
        'tumor_count', 'tumor_type_count', 'tumor_type_composition', 'indel_size',
        'transcript_id', 'transcript_id_version', 'type', 'version',
        'missense_count', 'trunc_count', 'inframe_count', 'splice_count',
        'total_count', 'missense_fraction', 'trunc_fraction', 'inframe_fraction',
        'splice_fraction',
    ]
    hotspots = hotspots[output_cols]

    remove_hotspots = (hotspots.type == "single residue") & ((hotspots.trunc_fraction > .75) | (hotspots.missense_count.isin([0,1])))
    n_removed = remove_hotspots.sum()
    if args.hotspots_2d_v3:
        # with v3 data, more hotspots may be filtered
        assert n_removed >= 120, f"Expected at least 120 removed hotspots, got {n_removed}"
    else:
        assert n_removed == 120, f"Expected 120 removed hotspots, got {n_removed}"
    print(f"Removing {n_removed} hotspots", file=sys.stderr)
    hotspots[~remove_hotspots].to_csv(sys.stdout, sep="\t", index=False)
    if args.removed_hotspots:
        hotspots[remove_hotspots].to_csv(args.removed_hotspots, sep="\t", index=False)