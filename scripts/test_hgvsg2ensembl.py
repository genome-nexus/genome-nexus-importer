#!/usr/bin/env python3
"""Tests for hgvsg2ensembl.

Run: python3 -m unittest discover -s scripts -p 'test_hgvsg2ensembl.py'

The cases marked REGRESSION encode behaviour where
``VEPController.hgvsgToEnsembl()`` in genome-nexus-vep is wrong; they are the
reason this converter exists rather than reusing that function.
"""

import io
import unittest

from hgvsg2ensembl import Deferred, Faidx, convert, minimise, run


class Substitution(unittest.TestCase):
    def test_snv(self):
        self.assertEqual(convert("1:g.100A>C"),
                         "1 100 100 N/C + 1:g.100A>C")

    def test_stated_ref_is_not_trusted(self):
        # --lookup_ref overwrites ref from the FASTA, so the asserted reference
        # base is deliberately replaced by the placeholder.
        self.assertTrue(convert("1:g.100A>C").split()[3].startswith("N/"))

    def test_real_grch38_egfr_l858r(self):
        self.assertEqual(convert("7:g.55191822T>G"),
                         "7 55191822 55191822 N/G + 7:g.55191822T>G")


class Deletion(unittest.TestCase):
    def test_bare_single_base(self):
        self.assertEqual(convert("1:g.100del"),
                         "1 100 100 N/- + 1:g.100del")

    def test_bare_range(self):
        self.assertEqual(convert("1:g.100_102del"),
                         "1 100 102 N/- + 1:g.100_102del")

    def test_stated_sequence_is_dropped(self):
        # REGRESSION: hgvsgToEnsembl emits "1 100 102 N/GGG +", which
        # --lookup_ref resolves to GGG/GGG -- a null variant, not a deletion.
        # Verified against gn-vep-111-grch38 that del / delGGG / delTTT all
        # produce byte-identical annotations, so dropping the sequence is
        # lossless.
        self.assertEqual(convert("1:g.100_102delGGG"),
                         "1 100 102 N/- + 1:g.100_102delGGG")

    def test_stated_sequence_dropped_single_base(self):
        self.assertEqual(convert("1:g.100delA"),
                         "1 100 100 N/- + 1:g.100delA")

    def test_identifier_keeps_the_original_key(self):
        # The 6th field becomes VEP's `id` and therefore the mongo `_id`, so it
        # must survive normalisation untouched.
        self.assertEqual(convert("1:g.100_102delGGG").split(" ", 5)[5],
                         "1:g.100_102delGGG")

    def test_sequence_length_may_contradict_the_range(self):
        # ~0.15% of delSEQ keys in production have a sequence whose length does
        # not match its coordinate range.  Dropping the sequence makes these
        # convert cleanly instead of handing VEP an inconsistent ref/range pair.
        self.assertEqual(convert("1:g.100_102delGG"),
                         "1 100 102 N/- + 1:g.100_102delGG")

    def test_long_real_world_deletion(self):
        key = "15:g.91014027_91015010del" + "G" * 984
        self.assertEqual(convert(key).split()[:5],
                         ["15", "91014027", "91015010", "N/-", "+"])


class DelIns(unittest.TestCase):
    def test_plain_delins(self):
        self.assertEqual(convert("1:g.100_102delinsAT"),
                         "1 100 102 N/AT + 1:g.100_102delinsAT")

    def test_single_base_delins(self):
        self.assertEqual(convert("1:g.100delinsAT"),
                         "1 100 100 N/AT + 1:g.100delinsAT")

    def test_compound_form_is_handled(self):
        # REGRESSION: hgvsgToEnsembl's trailing group is [A-Z]*, which cannot
        # match the lowercase "ins", so it rejects this outright -- 127,076 keys
        # in the production corpus.
        self.assertEqual(convert("1:g.100_102delGGGinsAT"),
                         "1 100 102 N/AT + 1:g.100_102delGGGinsAT")

    def test_inserted_sequence_is_never_dropped(self):
        for key in ("1:g.100_102delinsAT", "1:g.100_102delGGGinsAT"):
            self.assertEqual(convert(key).split()[3], "N/AT", key)


class Insertion(unittest.TestCase):
    def test_coordinates_are_swapped(self):
        # Ensembl format places an insertion by giving start = end + 1.
        self.assertEqual(convert("1:g.100_101insAT"),
                         "1 101 100 -/AT + 1:g.100_101insAT")

    def test_ref_is_dash_not_placeholder(self):
        self.assertEqual(convert("1:g.100_101insAT").split()[3], "-/AT")

    def test_inserted_sequence_is_never_dropped(self):
        self.assertEqual(convert("1:g.100_101insGGGG").split()[3], "-/GGGG")

    def test_non_adjacent_flanks_are_deferred(self):
        with self.assertRaises(Deferred):
            convert("1:g.100_105insAT")

    def test_single_position_insertion_is_deferred(self):
        # No second coordinate means the insertion point is ambiguous.
        with self.assertRaises(Deferred):
            convert("1:g.100insAT")


class NeedsReferenceSequence(unittest.TestCase):
    def test_dup_deferred(self):
        with self.assertRaises(Deferred):
            convert("1:g.100_102dup")

    def test_dup_with_sequence_deferred(self):
        with self.assertRaises(Deferred):
            convert("1:g.100_102dupGGG")

    def test_inv_deferred(self):
        with self.assertRaises(Deferred):
            convert("1:g.100_102inv")

    def test_deferred_not_dropped(self):
        # The whole point: hgvsgToEnsembl throws and the HTTP layer discards the
        # error when any sibling in the batch succeeded.  Here every rejection
        # must land in the deferral file.
        keys = ["1:g.100A>C", "1:g.100_102dup", "1:g.100_102inv", "1:g.100_102del"]
        out, defer = io.StringIO(), io.StringIO()
        stats, _ = run(io.StringIO("\n".join(keys)), out, defer)
        self.assertEqual(stats["read"], 4)
        self.assertEqual(stats["converted"], 2)
        self.assertEqual(stats["deferred"], 2)
        self.assertEqual(defer.getvalue().split(),
                         ["1:g.100_102dup", "1:g.100_102inv"])


class Malformed(unittest.TestCase):
    def test_no_colon(self):
        with self.assertRaises(Deferred):
            convert("BRAF")

    def test_not_genomic_notation(self):
        with self.assertRaises(Deferred):
            convert("ENST00000275493.7:c.2369C>T")

    def test_end_before_start(self):
        with self.assertRaises(Deferred):
            convert("1:g.102_100del")

    def test_empty_change(self):
        with self.assertRaises(Deferred):
            convert("1:g.100")

    def test_lowercase_sequence_is_deferred_not_guessed(self):
        with self.assertRaises(Deferred):
            convert("1:g.100a>c")


class Chromosomes(unittest.TestCase):
    def test_sex_and_mito(self):
        for chrom in ("X", "Y", "M", "MT"):
            self.assertEqual(convert("%s:g.100A>C" % chrom).split()[0], chrom)

    def test_assembly_contig_passes_through(self):
        # Not guessed at or rewritten; if VEP does not know it, the variant
        # becomes a straggler for the database pass.
        self.assertEqual(convert("KI270438.1:g.100A>C").split()[0], "KI270438.1")

    def test_dotted_contig_does_not_confuse_the_position_parser(self):
        self.assertEqual(convert("GL000205.2:g.4155del").split()[:4],
                         ["GL000205.2", "4155", "4155", "N/-"])


class Accounting(unittest.TestCase):
    def test_every_line_is_accounted_for(self):
        keys = ["1:g.100A>C", "1:g.100_102delGGG", "1:g.100_102dup",
                "", "  ", "1:g.100_101insAT", "garbage"]
        out, defer = io.StringIO(), io.StringIO()
        stats, _ = run(io.StringIO("\n".join(keys)), out, defer)
        self.assertEqual(stats["read"], 5)
        self.assertEqual(stats["converted"] + stats["deferred"], stats["read"])
        self.assertEqual(len(out.getvalue().strip().split("\n")), stats["converted"])

    def test_input_order_is_preserved(self):
        keys = ["3:g.300A>C", "1:g.100A>C", "2:g.200A>C"]
        out, defer = io.StringIO(), io.StringIO()
        run(io.StringIO("\n".join(keys)), out, defer)
        self.assertEqual([l.split()[0] for l in out.getvalue().strip().split("\n")],
                         ["3", "1", "2"])


if __name__ == "__main__":
    unittest.main(verbosity=2)


class Minimise(unittest.TestCase):
    """The trim itself, independent of any FASTA."""

    def test_shared_suffix(self):
        # 1:g.100011882_100011883delinsTT, ref CT -> database mode reports C/T
        self.assertEqual(minimise(100011882, 100011883, "CT", "TT"),
                         (100011882, 100011882, "C", "T"))

    def test_becomes_pure_deletion(self):
        # 1:g.100022945_100022946delinsA, ref GA -> G/-
        self.assertEqual(minimise(100022945, 100022946, "GA", "A"),
                         (100022945, 100022945, "G", "-"))

    def test_longer_shared_suffix(self):
        # 1:g.100015300_100015303delinsTTT, ref AGGT -> AGG/TT
        self.assertEqual(minimise(100015300, 100015303, "AGGT", "TTT"),
                         (100015300, 100015302, "AGG", "TT"))

    def test_shared_prefix(self):
        self.assertEqual(minimise(100, 102, "AGG", "AT"),
                         (101, 102, "GG", "T"))

    def test_becomes_insertion_start_exceeds_end(self):
        # Ensembl places an insertion with start = end + 1.
        s, e, r, a = minimise(100, 101, "AG", "AGT")
        self.assertEqual((s, e, r, a), (102, 101, "-", "T"))
        self.assertGreater(s, e)

    def test_nothing_to_trim(self):
        self.assertEqual(minimise(100, 102, "AGG", "TTC"),
                         (100, 102, "AGG", "TTC"))


class FaidxBackedDelins(unittest.TestCase):
    """End-to-end trimming through convert(), against a synthetic FASTA."""

    @classmethod
    def setUpClass(cls):
        import tempfile, os
        cls.dir = tempfile.mkdtemp()
        cls.fa = os.path.join(cls.dir, "ref.fa")
        # 60bp-per-line FASTA, chr1 = 120 bases of a known pattern
        seq = ("ACGTACGTAC" * 12)
        with open(cls.fa, "w") as fh:
            fh.write(">1\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
        with open(cls.fa + ".fai", "w") as fh:
            fh.write("1\t%d\t3\t60\t61\n" % len(seq))
        cls.seq = seq
        cls.fasta = Faidx(cls.fa)

    def test_faidx_reads_across_a_line_break(self):
        # bases 58-63 span the first newline
        self.assertEqual(self.fasta.fetch("1", 58, 63), self.seq[57:63])

    def test_faidx_single_base(self):
        self.assertEqual(self.fasta.fetch("1", 1, 1), self.seq[0])

    def test_faidx_out_of_range(self):
        self.assertIsNone(self.fasta.fetch("1", 119, 130))
        self.assertIsNone(self.fasta.fetch("nosuchchrom", 1, 2))

    def test_delins_trimmed_with_fasta(self):
        # ref at 1..2 is "AC"; delins "GC" shares the suffix C -> A/G at 1..1
        self.assertEqual(convert("1:g.1_2delinsGC", self.fasta),
                         "1 1 1 A/G + 1:g.1_2delinsGC")

    def test_delins_untrimmed_without_fasta(self):
        self.assertEqual(convert("1:g.1_2delinsGC"),
                         "1 1 2 N/GC + 1:g.1_2delinsGC")

    def test_compound_delins_trimmed(self):
        self.assertEqual(convert("1:g.1_2delACinsGC", self.fasta),
                         "1 1 1 A/G + 1:g.1_2delACinsGC")

    def test_other_classes_unaffected_by_fasta(self):
        # Only delins consults the FASTA; everything else is already minimal.
        for key in ("1:g.5A>C", "1:g.5_7del", "1:g.5_7delACG", "1:g.5_6insTT"):
            self.assertEqual(convert(key), convert(key, self.fasta), key)

    def test_unknown_chrom_falls_back_to_placeholder(self):
        self.assertEqual(convert("99:g.1_2delinsGC", self.fasta),
                         "99 1 2 N/GC + 99:g.1_2delinsGC")
