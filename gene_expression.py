import argparse
import polars as pl


def process_fetures(df):
    dfe = df.group_by(["gene","wellid"]).agg(pl.col("umi").n_unique())
    dfe = dfe.pivot(
    values="umi",
    index="gene",
    columns="wellid"
    )
    dfe = dfe.drop('null').filter(pl.col('gene') != 'null')
    dfe = dfe.fill_null(strategy='zero')
    return dfe


def main(args):
    df = pl.read_csv(
        source=args.input, 
        separator='\t', 
        has_header=True,
        columns=["gene", "wellid", "umi"]
    )
    df = process_fetures(df)
    df.write_csv(args.output, separator='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barcodes correction and assignment")
    parser.add_argument("-i", "--input", required=True, help="Input tsv file recording gene, cellid and umis")
    parser.add_argument("-o", "--output", required=True, help="Output count matrix in tsv")
    args = parser.parse_args()
    main(args)