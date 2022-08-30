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
import requests

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

    def test_new_transcript_id_is_valid(self):
        result = hotspots.update_hotspots_to_grch38.new_transcript_id_is_valid('ENST00000288602', 'ENST00000646891')
        self.assertTrue(result)

        result = hotspots.update_hotspots_to_grch38.new_transcript_id_is_valid('ENST00000275493', 'ENST00000646891')
        self.assertFalse(result)

    def test_get_translated_protein_sequence(self):
        result = hotspots.update_hotspots_to_grch38.get_translated_protein_sequence(hotspots.update_hotspots_to_grch38.ENSEMBL_GRCH37_SERVER,
                                                                                    'ENST00000288602')
        self.assertIsNotNone(result)
        result = hotspots.update_hotspots_to_grch38.get_translated_protein_sequence(hotspots.update_hotspots_to_grch38.ENSEMBL_GRCH38_SERVER,
                                                                                    'ENST_invalid')
        self.assertIsNone(result)
        
        # expect exception:
        with self.assertRaises(requests.exceptions.ConnectionError):
            hotspots.update_hotspots_to_grch38.get_translated_protein_sequence(hotspots.update_hotspots_to_grch38.ENSEMBL_GRCH38_SERVER+'err',
                                                                                    '')


    def test_generate_updated_grch38_hotspots_info(self):
        result_df = hotspots.update_hotspots_to_grch38.generate_updated_grch38_hotspots_info('data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt',
                                                                                             'mskcc',
                                                                                             'scripts/test_files/grch37_ensembl92/export/small_hotspots_v2_and_3d.txt')
                                                                                             #'data/grch37_ensembl92/export/hotspots_v2_and_3d.txt')
        # iterate over all rows and check if each expected change occured, e.g. for BRAF:
        for _, row in result_df.iterrows():
            hugo_symbol = row['hugo_symbol']
            if hugo_symbol == "BRAF":
                transcript_id = row['transcript_id']
                self.assertEqual('ENST00000646891', transcript_id)


    def test_find_grch38_transcript_id(self):
        hugo_and_transcript_map = hotspots.update_hotspots_to_grch38.get_gene_and_transcript_map('data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt')
        # or test wit a mock:
        #  hugo_and_transcript_map = {
        #     'BRAF': {'previous_symbols': 'BLA,ETC', 'synonyms': None},
        #     'NEW_H3F3A': {'previous_symbols': 'BLA, H3F3A', 'synonyms': None, 'mskcc_canonical_transcript': 'ENST00000366815'},
        # }
        transcript_id = hotspots.update_hotspots_to_grch38.find_grch38_transcript_id('H3F3A', hugo_and_transcript_map, 'mskcc')
        self.assertEqual('ENST00000366815', transcript_id)

        transcript_id = hotspots.update_hotspots_to_grch38.find_grch38_transcript_id('invalidhugo', hugo_and_transcript_map, 'mskcc')
        self.assertIsNone(transcript_id)
