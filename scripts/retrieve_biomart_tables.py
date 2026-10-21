#!/usr/bin/env python3
"""Default way to retrieve BioMart tables for homo_sapiens / GRCh37 or GRCh38.

Connects directly to Ensembl's public, anonymous, read-only MySQL mirror
(ensembldb.ensembl.org), which serves every release's core database
independently of the Ensembl website/BioMart REST layer (which goes down
periodically, especially for archived historical releases), and reproduces
the same 4 BioMart input files (geneids, pfam, refseq, ccds) that
retrieve_biomart_tables.R would otherwise fetch via biomaRt.

Usage: retrieve_biomart_tables.py <output_dir> <genome_build> <ensembl_release>
Example: retrieve_biomart_tables.py grch38_ensembl111/input/ grch38 111

Note: GRCh37 was frozen at Ensembl release 75 (homo_sapiens_core_75_37) -
there is no core database for later release numbers under GRCh37, so
<ensembl_release> must be 75 (or earlier) when <genome_build> is grch37.
"""

import subprocess
import sys

MYSQL_HOST = "ensembldb.ensembl.org"
MYSQL_PORT = "3306"
MYSQL_USER = "anonymous"

ASSEMBLY_SUFFIX = {"grch37": "37", "grch38": "38"}

QUERIES = {
    "ensembl_biomart_geneids.txt": (
        ["Gene stable ID", "Transcript stable ID", "Versioned transcript ID", "HGNC symbol", "HGNC ID"],
        """
        SELECT g.stable_id, t.stable_id, CONCAT(t.stable_id, '.', t.version),
               IFNULL(gh.display_label, ''), IFNULL(REPLACE(gh.dbprimary_acc, 'HGNC:', ''), '')
        FROM gene g
        JOIN transcript t ON t.gene_id = g.gene_id
        LEFT JOIN (
            SELECT ox.ensembl_id AS gene_id, x.display_label, x.dbprimary_acc
            FROM object_xref ox
            JOIN xref x ON x.xref_id = ox.xref_id
            JOIN external_db ed ON ed.external_db_id = x.external_db_id AND ed.db_name = 'HGNC'
            WHERE ox.ensembl_object_type = 'Gene'
        ) gh ON gh.gene_id = g.gene_id
        WHERE g.is_current = 1 AND t.is_current = 1;
        """,
    ),
    "ensembl_biomart_pfam.txt": (
        ["Gene stable ID", "Transcript stable ID", "Versioned transcript ID", "Gene name",
         "Pfam domain ID", "Pfam domain start", "Pfam domain end"],
        """
        SELECT g.stable_id, t.stable_id, CONCAT(t.stable_id, '.', t.version), IFNULL(gn.display_label, ''),
               IFNULL(pf.hit_name, ''), IFNULL(pf.seq_start, ''), IFNULL(pf.seq_end, '')
        FROM gene g
        JOIN transcript t ON t.gene_id = g.gene_id
        JOIN translation tl ON tl.transcript_id = t.transcript_id
        LEFT JOIN xref gn ON gn.xref_id = g.display_xref_id
        LEFT JOIN (
            SELECT pf.translation_id, pf.hit_name, pf.seq_start, pf.seq_end
            FROM protein_feature pf
            JOIN analysis a ON a.analysis_id = pf.analysis_id AND a.logic_name = 'pfam'
        ) pf ON pf.translation_id = tl.translation_id
        WHERE g.is_current = 1 AND t.is_current = 1;
        """,
    ),
    "ensembl_biomart_refseq.txt": (
        ["Transcript stable ID", "Versioned transcript ID", "RefSeq mRNA ID"],
        """
        SELECT t.stable_id, CONCAT(t.stable_id, '.', t.version), IFNULL(rs.dbprimary_acc, '')
        FROM transcript t
        LEFT JOIN (
            SELECT ox.ensembl_id AS transcript_id, x.dbprimary_acc
            FROM object_xref ox
            JOIN xref x ON x.xref_id = ox.xref_id
            JOIN external_db ed ON ed.external_db_id = x.external_db_id
                AND ed.db_name IN ('RefSeq_mRNA', 'RefSeq_mRNA_predicted')
            WHERE ox.ensembl_object_type = 'Transcript'
        ) rs ON rs.transcript_id = t.transcript_id
        WHERE t.is_current = 1;
        """,
    ),
    "ensembl_biomart_ccds.txt": (
        ["Transcript stable ID", "Versioned transcript ID", "CCDS ID"],
        """
        SELECT t.stable_id, CONCAT(t.stable_id, '.', t.version), IFNULL(cc.dbprimary_acc, '')
        FROM transcript t
        LEFT JOIN (
            SELECT ox.ensembl_id AS transcript_id, x.dbprimary_acc
            FROM object_xref ox
            JOIN xref x ON x.xref_id = ox.xref_id
            JOIN external_db ed ON ed.external_db_id = x.external_db_id AND ed.db_name = 'CCDS'
            WHERE ox.ensembl_object_type = 'Transcript'
        ) cc ON cc.transcript_id = t.transcript_id
        WHERE t.is_current = 1;
        """,
    ),
}


def run_query(database, sql):
    result = subprocess.run(
        ["mysql", "-h", MYSQL_HOST, "-P", MYSQL_PORT, "-u", MYSQL_USER, "-N", "-B", database, "-e", sql],
        capture_output=True, text=True, check=True,
    )
    return result.stdout


def database_exists(database):
    result = subprocess.run(
        ["mysql", "-h", MYSQL_HOST, "-P", MYSQL_PORT, "-u", MYSQL_USER, "-N", "-B",
         "-e", f"SHOW DATABASES LIKE '{database}';"],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip() == database


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <output_dir> <genome_build> <ensembl_release>", file=sys.stderr)
        print(f"Example: {sys.argv[0]} grch38_ensembl111/input/ grch38 111", file=sys.stderr)
        sys.exit(1)

    output_dir, genome_build, release = sys.argv[1], sys.argv[2], sys.argv[3]

    if genome_build not in ASSEMBLY_SUFFIX:
        print(f"Unsupported genome_build '{genome_build}', expected one of {sorted(ASSEMBLY_SUFFIX)}", file=sys.stderr)
        sys.exit(1)

    database = f"homo_sapiens_core_{release}_{ASSEMBLY_SUFFIX[genome_build]}"

    if not database_exists(database):
        print(f"No core database '{database}' on {MYSQL_HOST} for {genome_build} release {release}.", file=sys.stderr)
        if genome_build == "grch37":
            print("GRCh37 was frozen at Ensembl release 75 - use release 75 or earlier.", file=sys.stderr)
        sys.exit(1)

    for filename, (header, sql) in QUERIES.items():
        print(f"Querying {database} for {filename}...")
        tsv_rows = run_query(database, sql)
        out_path = f"{output_dir.rstrip('/')}/{filename}"
        with open(out_path, "w") as f:
            f.write("\t".join(header) + "\n")
            f.write(tsv_rows)
        print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()
