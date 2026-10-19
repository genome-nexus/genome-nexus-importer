#!/usr/bin/env python2.7

"""
Copyright (c) 2018 The Hyve B.V.
This code is licensed under the GNU Affero General Public License (AGPL),
version 3, or (at your option) any later version.
"""

import unittest
import difflib
import os
import transform_gff_to_tsv_for_exon_info_from_ensembl
import gzip
import hotspots.hotspot_transcript_utils

class TransformTestCase(unittest.TestCase):
    """Superclass for testcases that test the transformation steps.
    """

    def assertFileGenerated(self, tmp_file_name, expected_file_name):
        """Assert that a file has been generated with the expected contents."""
        self.assertTrue(os.path.exists(tmp_file_name))
        with gzip.open(tmp_file_name, 'rU') as out_file, gzip.open(expected_file_name, 'rU') as ref_file:
            base_filename = os.path.basename(tmp_file_name)
            base_input = os.path.basename(expected_file_name)
            diff_result = difflib.context_diff(
                ref_file.readlines(),
                out_file.readlines(),
                fromfile='Expected {}'.format(base_input),
                tofile='Generated {}'.format(base_filename))
        diff_line_list = list(diff_result)
        self.assertEqual(diff_line_list, [], msg='\n' + ''.join(diff_line_list))
        # remove temp file if all is fine:
        try:
            os.remove(tmp_file_name)
        except WindowsError:
            # ignore this Windows specific error...probably happens because of virus scanners scanning the temp file...
            pass

    def test_pfam_transformation_step(self):
        """Test pfam TSV to internal data structure transformation"""
        # TODO

    def test_gff_to_tsv(self):
        """Test transcript info gff to internal data structure transformation"""
        # Build up arguments and run
        out_file_name = 'test_files/transform_gff_to_tsv/ensembl_transcript_info.txt.gz~'
        gff_input_file_name = 'test_files/transform_gff_to_tsv/sub_Homo_Sapiens.gff3.gz'
        transform_gff_to_tsv_for_exon_info_from_ensembl.transform_gff_to_tsv(gff_input_file_name, out_file_name)
        self.assertFileGenerated(out_file_name, 'test_files/transform_gff_to_tsv/ensembl_transcript_info.txt.gz')

    # The following tests cover hotspot_transcript_utils.py, the alignment/translation
    # logic shared by overwrite_transcript_and_residue_to_mskcc_isoform.py (GRCh37 native
    # transcript -> MSKCC's preferred GRCh37 transcript) and liftover_hotspots_to_grch38.py
    # (GRCh37 -> GRCh38). These replace the old tests for the now-removed
    # update_hotspots_to_grch38.py, which depended on live Ensembl REST API calls
    # (grch37.rest.ensembl.org / rest.ensembl.org) rather than release-pinned sequence data,
    # and dropped a gene's hotspots entirely on any whole-protein sequence mismatch instead
    # of translating per-residue.

    def test_build_position_map_identical_sequences(self):
        seq = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq, seq)
        self.assertEqual(1, pos_map[1])
        self.assertEqual(len(seq), pos_map[len(seq)])

    def test_build_position_map_with_internal_deletion(self):
        # seq_b is seq_a with a 5-residue internal deletion (mirrors the real RAC1 case:
        # two transcripts identical except for one internal indel, so positions before the
        # deletion are unchanged and positions after it shift by a constant amount).
        seq_a = 'MQAIKXXXXXCVVVGDGAVGKT'
        seq_b = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq_a, seq_b)
        self.assertEqual(5, pos_map[5])  # before the deletion: unchanged
        self.assertIsNone(pos_map[7])  # inside the deleted region: no corresponding position
        self.assertEqual(10, pos_map[15])  # after the deletion: shifted by -5

    def test_translate_residue_single_residue_unchanged(self):
        seq = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq, seq)
        fields, status, detail = hotspots.hotspot_transcript_utils.translate_residue(
            'M1', 'single residue', pos_map, seq, seq)
        self.assertEqual({'residue': 'M1', 'amino_acid_position': '1'}, fields)
        self.assertEqual('unchanged', status)

    def test_translate_residue_single_residue_shifted(self):
        seq_a = 'MQAIKXXXXXCVVVGDGAVGKT'
        seq_b = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq_a, seq_b)
        # position 15 in seq_a ('G') is the same residue as position 10 in seq_b
        fields, status, detail = hotspots.hotspot_transcript_utils.translate_residue(
            'G15', '3d', pos_map, seq_a, seq_b)
        self.assertEqual({'residue': 'G10', 'amino_acid_position': '10'}, fields)
        self.assertEqual('shifted', status)

    def test_translate_residue_deleted(self):
        seq_a = 'MQAIKXXXXXCVVVGDGAVGKT'
        seq_b = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq_a, seq_b)
        # position 7 falls inside the region deleted in seq_b
        fields, status, detail = hotspots.hotspot_transcript_utils.translate_residue(
            'X7', 'single residue', pos_map, seq_a, seq_b)
        self.assertIsNone(fields)
        self.assertEqual('deleted', status)

    def test_translate_residue_ref_mismatch(self):
        seq_a = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq_a, seq_a)
        # claims position 1 is 'Z', but it's actually 'M' - the hotspot's own stated
        # reference amino acid doesn't even match the source sequence
        fields, status, detail = hotspots.hotspot_transcript_utils.translate_residue(
            'Z1', 'single residue', pos_map, seq_a, seq_a)
        self.assertIsNone(fields)
        self.assertEqual('grch37_ref_mismatch', status)

    def test_translate_residue_in_frame_indel_range(self):
        seq_a = 'MQAIKXXXXXCVVVGDGAVGKT'
        seq_b = 'MQAIKCVVVGDGAVGKT'
        pos_map = hotspots.hotspot_transcript_utils.build_position_map(seq_a, seq_b)
        fields, status, detail = hotspots.hotspot_transcript_utils.translate_residue(
            '1-3', 'in-frame indel', pos_map, seq_a, seq_b)
        self.assertEqual({'residue': '1-3', 'amino_acid_position': '1-3'}, fields)
        self.assertEqual('unchanged', status)
