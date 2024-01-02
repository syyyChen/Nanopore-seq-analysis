import argparse
import os
import polars as pl


def process_fetures(tsv):
    dfe = pl.read_csv(
        source=tsv, 
        separator='\t', 
        has_header=True,
        columns=["gene", "wellid", "umi"]
    )
    dfe = dfe.group_by(["gene","wellid"]).agg(pl.col("umi").n_unique())
    dfe = dfe.pivot(
    values="umi",
    index="gene",
    columns="wellid"
    )
    dfe = dfe.drop('null').filter(pl.col('gene') != 'null')
    dfe = dfe.fill_null(strategy='zero')
    return dfe


def create_meta(tsv):
    df = pl.read_csv(
        source=tsv, 
        separator='\t', 
        has_header=True,
        columns=["assign_status", "wellid", "pattern", "umi"]
    )
    meta = df.group_by('wellid').agg([
            pl.count().alias("total_reads"),
            pl.col("assign_status").filter(pl.col("assign_status") != "Unassigned_Unmapped").count().alias("mapped"),
            pl.col("assign_status").filter(pl.col("assign_status") == "Assigned").count().alias("assigned"),
            pl.col("pattern").filter(pl.col("pattern") == "+").count().alias("positive_pattern"),
            pl.col("umi").filter(pl.col("umi").is_not_null()).count().alias("qualified_umi")
            ])
    meta = meta.with_columns([
            (pl.col("mapped") / pl.col("total_reads") * 100).alias("pct_mapped_reads"),
            (pl.col("assigned") / pl.col("total_reads") * 100).alias("pct_assigned_reads"),
            (pl.col("qualified_umi") / pl.col("total_reads") * 100).alias("pct_qualified_umi"),
            (pl.col("positive_pattern") / pl.col("total_reads") * 100).alias("pct_positive_pattern")
        ]).select(
            ["wellid", "total_reads", "pct_mapped_reads", "pct_assigned_reads", "pct_qualified_umi", "pct_positive_pattern"]
            ).filter(pl.col('wellid') != 'null')
    return meta


def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    count_matrix = os.path.join(args.output_dir, "gene_count.tsv")
    metadata = os.path.join(args.output_dir, "meta.tsv")
    
    ge = process_fetures(args.input)
    meta = create_meta(args.input)
    
    ge.write_csv(count_matrix, separator='\t')
    meta.write_csv(metadata, separator='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barcodes correction and assignment")
    parser.add_argument("-i", "--input", required=True, help="Input tsv file recording gene, cellid, assigned status, pattern and umis")
    parser.add_argument("-o", "--output_dir", required=True, help="Output path")
    args = parser.parse_args()
    main(args)