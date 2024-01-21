import polars as pl
import argparse


def ascii_to_phred(ascii_str):
    return [ord(char) - 33 for char in ascii_str]

def process_umis(umis, qual_threshold):
    umi_df = pl.read_csv(
        source=umis, 
        separator='\t', 
        has_header=True
        )
    umi_df = umi_df.with_columns(
        pl.col("UR_qual")
        .map_elements(ascii_to_phred)
        .alias("phred_scores")
        )
    umi_df = umi_df.with_columns(
        avg_phred_score=pl.col("phred_scores").list.mean()
        )
    umi_df = umi_df.with_columns(
        (pl.col("avg_phred_score") >= qual_threshold)
        .alias("qualified")
        )
    umi_df = umi_df.with_columns(
        pl.when(pl.col("qualified"))
        .then(pl.col("UR"))
        .otherwise(pl.lit(None))
        .alias("UR")
    )
    umi_df = umi_df.select(pl.col(["read_id", "UR"]))
    return umi_df


def main(args):
    
    gene_df = pl.read_csv(source=args.gene_assign, 
            separator='\t', 
            has_header=False,
            new_columns=["read_id", "assign_status", "assign_n", "gene"])
    gene_df = gene_df.with_columns(pl.col("gene").map_elements(lambda x: None if x == "NA" else x))
    
    bc_df = pl.read_csv(source=args.barcodes, 
            separator='\t', 
            has_header=True,
            columns=["read_id", "wellid", "pattern"])
    bc_df = bc_df.with_columns(pl.col("wellid").map_elements(lambda x: None if x == "Unsure" else x))
    
    umi_df = process_umis(
        umis=args.umis, 
        qual_threshold=args.qual_threshold)
    
    df = umi_df.join(gene_df, on="read_id", how="inner").join(bc_df, on="read_id", how="inner")
    
    df.write_csv(file=args.output,separator='\t')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge tags for UMI correction")
    parser.add_argument("-o", "--output", required=True, help="Output tsv file")
    parser.add_argument("-g", "--gene_assign", required=True, help="Gene assignment file from featureCounts")
    parser.add_argument("-b", "--barcodes", required=True, help="Corrected barcodes and cell ids")
    parser.add_argument("-u", "--umis", required=True, help="Uncorrected UMIs")
    parser.add_argument("-q", "--qual_threshold", type=int, default=10, help="Minimum average quality score for UMIs filtering")
    args = parser.parse_args()
    main(args)