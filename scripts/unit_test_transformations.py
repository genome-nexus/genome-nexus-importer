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
import hotspots.update_hotspots_to_grch38


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

    def test_get_gene_and_transcript_map(self):
        hugo_and_transcript_map = hotspots.update_hotspots_to_grch38.get_gene_and_transcript_map('data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt')
        self.assertEqual('ENST00000646891', hugo_and_transcript_map['BRAF']['mskcc_canonical_transcript'])
        self.assertEqual('ENST00000288602', hugo_and_transcript_map['BRAF']['uniprot_canonical_transcript'])

    def test_validate_new_transcript_id(self):
        result = hotspots.update_hotspots_to_grch38.validate_new_transcript_id('ENST00000288602', 'ENST00000646891')
        self.assertTrue(result)

        result = hotspots.update_hotspots_to_grch38.validate_new_transcript_id('ENST00000275493', 'ENST00000646891')
        self.assertFalse(result)
