import os
import argparse
import sys
import polars as pl
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

vsearch_colnames = [
    "query",
    "target",
    "id",
    "alnlen",
    "mism",
    "opens",
    "qilo",
    "qihi",
    "qstrand",
    "tilo",
    "tihi",
    "ql",
    "tl",
]

schema = {
    'query': pl.Utf8,
    'target': pl.Utf8,
    'id': pl.Float32,
    'alnlen': pl.UInt32,
    'mism': pl.UInt32,
    'opens': pl.UInt32,
    'qilo': pl.UInt32,
    'qihi': pl.UInt32,
    'qstrand': pl.Utf8,
    'tilo': pl.UInt32,
    'tihi': pl.UInt32,
    'ql': pl.UInt32,
    'tl': pl.UInt32
}

valid_patterns = (("adapter1_f", "adapter2_f"), ("adapter2_r", "adapter1_r"))



def read_vsearch_result(tsv_file, selected_columns=vsearch_colnames):
    """Read vsearch results as a polars dataframe"""
    for col in selected_columns:
        if col not in schema:
            raise ValueError(f"Invalid column name: {col}")
    
    columns_idxs = [vsearch_colnames.index(x) for x in selected_columns]
    selected_schema = {col: schema[col] for col in selected_columns}
    
    df = pl.read_csv(
        source=tsv_file,
        has_header=False,
        separator='\t',
        dtypes=selected_schema,
        columns=columns_idxs,
        new_columns=list(selected_schema.keys())
    )
    return df



def find_valid_pairs(vs_result):
    print("Finding valid pairs...")
    
    # Create shifted columns for targets and positions
    vs_result = vs_result.select(
        ['query','target','qilo','qihi']
        ).with_columns([
        vs_result['target'].shift(-1).alias('next_target'),
        vs_result['qihi'].shift(-1).alias('next_qihi'),
        vs_result['qilo'].shift(-1).alias('next_qilo')
    ])
    
    vs_result = vs_result.with_columns(
        pl.when((pl.col("target") == valid_patterns[0][0]) & (pl.col("next_target") == valid_patterns[0][1])).then(pl.lit("+"))
        .when((pl.col("target") == valid_patterns[1][0]) & (pl.col("next_target") == valid_patterns[1][1])).then(pl.lit("-"))
        .otherwise(pl.lit("*")).alias("pattern")
    )
    
    # Filter for valid pairs
    valid_pairs_condition = ((vs_result['pattern'] != '*') & (vs_result['qihi'] < vs_result['next_qilo'] - 10))
    valid_pairs_df = vs_result.filter(valid_pairs_condition)
    
    valid_pairs_df = valid_pairs_df.with_columns([
        pl.when(pl.col('pattern') == '+')
        .then(pl.col('qilo') - 11)
        .otherwise(pl.col('next_qihi'))
        .alias('i7_start'),

        pl.when(pl.col('pattern') == '+')
        .then(pl.col('next_qihi'))
        .otherwise(pl.col('qilo') - 11)
        .alias('i5_start')
        ])
    
    # Group and aggregate
    valid_query = valid_pairs_df.group_by('query').agg([
        pl.col('qihi').alias('start_positions'),
        pl.col('next_qilo').alias('end_positions'),
        pl.col('pattern').alias('patterns'),
        pl.col('i7_start').alias('i7_start'),
        pl.col('i5_start').alias('i5_start')
    ])
    valid_query = valid_query.with_columns(pl.col("patterns").list.len().alias("n_reads"))

    # Create dictionary of query to valid pairs
    valid_pairs_dict = {}
    for row in valid_query.iter_rows():
        query_id = row[0]
        start_positions = row[1]
        end_positions = [x - 1 for x in row[2]]
        patterns = row[3]
        i7 = row[4]
        i5 = row[5]
        valid_pairs_dict[query_id] = list(zip(start_positions, end_positions, patterns, i7, i5))
    
    return valid_pairs_dict



def write_barcode_batch(bc_batch, bc_output, write_header=False):
    with open(bc_output, 'a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["read_id", "pattern", "i7index_uncorr", "i5index_uncorr"], delimiter='\t')
        if write_header:
            writer.writeheader()
        writer.writerows(bc_batch)


def extract_barcodes(input_bam, valid_pairs_dict, output_dir, flank_in = 5, flank_out = 5, batch_size=1000):
    """Extract reads from the bam file with batch writing."""
    print("Extracting barcodes...")
    n_reads = 0
    i7_output = os.path.join(output_dir, "index_i7.fastq")
    i5_output = os.path.join(output_dir, "index_i5.fastq")
    bc_tsv = os.path.join(output_dir, "barcodes.tsv")

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bamfile, \
        open(i7_output, 'a') as i7_file, open(i5_output, 'a') as i5_file:
            i7_batch = []
            i5_batch = []
            bc_batch = []
            
            for read in bamfile.fetch(until_eof=True):
                query_name = read.query_name
                if query_name in valid_pairs_dict:
                    multi_reads = 1 if len(valid_pairs_dict[query_name])> 1 else 0
                    for idx, (start, end, pattern, i7_start, i5_start) in enumerate(valid_pairs_dict[query_name]):
                        query_id = f"{query_name}_{idx}" if multi_reads else query_name
                        
                        # output the fastq file for some analysis
                        i7_end = min(i7_start + 10 + flank_in, read.rlen)
                        i7_start = max(i7_start - flank_out, 0)
                        i5_end = min(i5_start + 10 + flank_in, read.rlen)
                        i5_start = max(i5_start - flank_out, 0)
                        i7_seq = read.seq[i7_start:i7_end] 
                        i7_qual = read.qual[i7_start:i7_end]
                        i5_seq = read.seq[i5_start:i5_end]
                        i5_qual = read.qual[i5_start:i5_end]
                        i7_record = SeqRecord(Seq(i7_seq), id=query_id, description=pattern, \
                            letter_annotations={"phred_quality": [ord(q) - 33 for q in i7_qual]})
                        i5_record = SeqRecord(Seq(i5_seq), id=query_id, description=pattern, \
                            letter_annotations={"phred_quality": [ord(q) - 33 for q in i5_qual]})
                        i7_batch.append(i7_record)
                        i5_batch.append(i5_record)
                        
                        # output barcodes table for correction
                        barcode_info = {
                            "read_id": query_id,
                            "pattern": pattern,
                            "i7index_uncorr": i7_seq,
                            "i5index_uncorr": i5_seq
                            }
                        bc_batch.append(barcode_info)

                        n_reads += 1
                        if n_reads >= batch_size:
                            SeqIO.write(i7_batch, i7_file, "fastq")
                            SeqIO.write(i5_batch, i5_file, "fastq")
                            i7_batch.clear()
                            i5_batch.clear()
                            write_barcode_batch(bc_batch, bc_tsv, write_header=(n_reads == batch_size))
                            bc_batch.clear()
                            
            if bc_batch:
                SeqIO.write(i7_batch, i7_file, "fastq")
                SeqIO.write(i5_batch, i5_file, "fastq")
                write_barcode_batch(bc_batch, bc_tsv)



def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    vsearch_results = read_vsearch_result(args.vsearch_result)

    valid_pairs_dict = find_valid_pairs(vsearch_results)
    
    extract_barcodes(
        input_bam=args.input_bam,
        valid_pairs_dict=valid_pairs_dict,
        output_dir=args.output_dir,
        flank_in=args.flank_in,
        flank_out=args.flank_out
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads from basecalled bam file and collect statistics.")
    parser.add_argument("-b", "--input_bam", required=True, help="Input BAM file path")
    parser.add_argument("-v", "--vsearch_result", required=True, help="Input vsearch_result file path")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-fi", "--flank_in", type=int, default=5, help="Extra nucleotides to extract at two ends of the reads")
    parser.add_argument("-fo", "--flank_out", type=int, default=5, help="Extra nucleotides to extract stretching beyond the assumed index boundaries")
    args = parser.parse_args()
    main(args)