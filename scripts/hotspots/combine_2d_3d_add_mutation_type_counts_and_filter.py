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
    parser.add_argument("hotspots_2d", default="../data/hotspots/v2_multi_type_residue.txt", type=str, help="2D cancerhotspots data file")
    parser.add_argument("hotspots_3d", default="../data/hotspots/3d_hotspots.txt", type=str, help="3D cancerhotspots data file")
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

    hotspots_3d = pd.read_csv(args.hotspots_3d, sep="\t", dtype=str)
    hotspots_3d.columns = [c.lower().replace("-","_") for c in hotspots_3d.columns]
    # add type column
    hotspots_3d['type'] = '3d'

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
    
    remove_hotspots = (hotspots.type == "single residue") & ((hotspots.trunc_fraction > .75) | (hotspots.missense_count.isin([0,1])))
    # removing 120 hotspots
    assert(remove_hotspots.sum() == 120)
    hotspots[~remove_hotspots].to_csv(sys.stdout, sep="\t", index=False)
    if args.removed_hotspots:
        hotspots[remove_hotspots].to_csv(args.removed_hotspots, sep="\t", index=False)