import os
import argparse
import sys
import csv
import numpy as np
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from collections import defaultdict
import tempfile
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed


complement_trans = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)

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

def write_adapters_fasta(adapter1_seq, adapter2_seq, output):
    """Write adapters fasta for use with vsearch."""
    adapters = []
    for adapter, seq in {
        "adapter1_f": adapter1_seq,
        "adapter1_r": adapter1_seq[::-1].translate(complement_trans),
        "adapter2_f": adapter2_seq,
        "adapter2_r": adapter2_seq[::-1].translate(complement_trans),
    }.items():
        entry = SeqRecord(Seq(seq), id=adapter, name="", description="")

        adapters.append(entry)
        with open(output, 'w') as fh:
            SeqIO.write(adapters, fh, "fasta")


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


def mask_regions(vs_result_df, fasta, tmp_dir, batch_size=1000):
    """Mask mapped regions of reads with 'U' according to vsearch results."""
    regions_to_mask = defaultdict(list)
    for row in vs_result_df.iter_rows():
        if row[11] == 0:
            continue
        query_id = row[0]
        start = row[6] - 1
        end = row[7]
        regions_to_mask[query_id].append((start, end))

    masked_records = []
    with open(os.path.join(tmp_dir, "masked.fasta"), 'w') as out_fh:
        for record in SeqIO.parse(fasta, "fasta"):
            masked_seq_str = str(record.seq)
            for start, end in regions_to_mask[record.id]:
                masked_seq_str = masked_seq_str[:start] + 'U' * (end - start) + masked_seq_str[end:]
            masked_record = SeqRecord(Seq(masked_seq_str), id=record.id, description=record.description)
            masked_records.append(masked_record)

            if len(masked_records) >= batch_size:
                SeqIO.write(masked_records, out_fh, "fasta")
                masked_records.clear()

        if masked_records:
            SeqIO.write(masked_records, out_fh, "fasta")



def iterative_vsearch(input_bam, adapters_fasta, output_dir, tmp_dir, min_len=20, id=0.7, rounds=5, threads=1):
    """Do iterative vsearch to find all adapters in the reads for pattern identification"""
    
    output = os.path.join(output_dir, "adapter_scan.vsearch.tsv")
    res_all = pl.DataFrame(schema=schema)
    query = os.path.join(tmp_dir, "input.fasta")
    subprocess.run(f"samtools fasta {input_bam} > {query}", shell=True)
        
    for round in range(1, rounds+1):
        result = os.path.join(tmp_dir, f"vsearch_res_{round}.tsv")
        tmp_fasta = os.path.join(tmp_dir, "tmp.fasta")
        vsearch_cmd = f"vsearch --usearch_global {query} --db {adapters_fasta} \
        --minseqlength {min_len} --maxaccepts 5 --id {id} --strand plus \
        --wordlength 3 --minwordmatches 10 --threads {threads} --userfields \
        'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
        --userout {result} --matched {tmp_fasta}"
        if round == 1:
            vsearch_cmd += " --output_no_hits"
            
        subprocess.run(vsearch_cmd, shell=True, check=True)
        
        if not os.path.getsize(result):
            break
        vsearch_results = read_vsearch_result(result)
        res_all = pl.concat([res_all, vsearch_results])
        if round < rounds:
            mask_regions(vsearch_results, tmp_fasta, tmp_dir)
            query = os.path.join(tmp_dir,"masked.fasta")
    
    res_all = res_all.sort(["query", "qilo"])
    res_all.write_csv(output, separator="\t", has_header=False)
    return res_all


def find_valid_pairs(vs_result):
    print("Finding valid pairs...")
    stat = {}
    stat['no_adapter'] = vs_result.filter(pl.col('target') == '*').height
    n_with_adapter = vs_result.filter(pl.col('target') != '*').n_unique(subset = ['query'])
    
    # Create shifted columns for targets and positions
    vs_result = vs_result.select(
        ['query','target','qilo','qihi']
        ).with_columns([
        vs_result['query'].shift(-1).alias('next_query'),
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
    # Note: do not do subtraction with qilo & qihi here!
    valid_pairs_condition = ((vs_result['query'] == vs_result['next_query']) & 
                            (vs_result['pattern'] != '*') & 
                            (vs_result['next_qilo'] > vs_result['qihi'] + 30))
    valid_pairs_df = vs_result.filter(valid_pairs_condition)
    
    # Group and aggregate
    valid_query = valid_pairs_df.group_by('query').agg([
        pl.col('qihi').alias('start_positions'),
        pl.col('next_qilo').alias('end_positions'),
        pl.col('pattern').alias('patterns')
    ])
    
    valid_query = valid_query.with_columns(pl.col("patterns").list.len().alias("n_reads"))

    n_valid = valid_query.height
    stat['invalid'] = n_with_adapter - n_valid
    stat['one_read'] = valid_query.filter(pl.col("n_reads") == 1).height
    stat['two_reads'] = valid_query.filter(pl.col("n_reads") == 2).height
    stat['more_reads'] = n_valid - stat['one_read'] - stat['two_reads']

    # Create dictionary of query to valid pairs
    valid_pairs_dict = {}
    for row in valid_query.iter_rows():
        query_id = row[0]
        start_positions = row[1]
        end_positions = [x - 1 for x in row[2]]
        patterns = row[3]
        valid_pairs_dict[query_id] = list(zip(start_positions, end_positions, patterns))
    
    print(f"Found {n_valid} ONT reads with valid pairs.")
    return valid_pairs_dict, stat



def process_read(read, valid_pairs_dict):
    """Process a single Nanopore read for extraction."""
    query_name = read.query_name
    if query_name not in valid_pairs_dict:
        return None

    multi_reads = 1 if len(valid_pairs_dict[query_name]) > 1 else 0
    extracted_reads = []

    for idx, (start, end, pattern) in enumerate(valid_pairs_dict[query_name]):
        query_id = f"{query_name}_{idx}" if multi_reads else query_name
        new_read = pysam.AlignedSegment()
        new_read.query_name = query_id
        new_read.flag = read.flag
        seq = read.query_sequence[start:end]
        qual = read.query_qualities[start:end]
        new_read.query_sequence = seq if pattern == '+' else seq[::-1].translate(complement_trans)
        new_read.query_qualities = qual if pattern == '+' else qual[::-1]
        new_read.tags = read.tags 
        new_read.set_tag('mp', multi_reads, 'i')
        new_read.set_tag('pn', pattern, 'Z')
        extracted_reads.append(new_read)

    return extracted_reads



def concurrent_extract_reads(input_bam, valid_pairs_dict, output_bam, batch_size=1000, threads=1):
    """Extract reads concurrently and write to a BAM file."""
    n_reads = 0
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bamfile, \
        pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as outfile, \
        ThreadPoolExecutor(max_workers=threads) as executor:

        future_to_read = {
            executor.submit(process_read, read, valid_pairs_dict): 
                read for read in bamfile.fetch(until_eof=True)
                }
        reads_batch = []

        for future in as_completed(future_to_read):
            extracted_reads = future.result()
            if extracted_reads:
                reads_batch.extend(extracted_reads)
                n_reads += len(extracted_reads)

                if len(reads_batch) >= batch_size:
                    for b_read in reads_batch:
                        outfile.write(b_read)
                    reads_batch.clear()
                    

        # Write remaining reads in the last batch
        for b_read in reads_batch:
            outfile.write(b_read)
            
    print(f"{n_reads} reads with pattern have been extracted.")



def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    with tempfile.TemporaryDirectory() as tmp_dir:
        adapters_fasta = os.path.join(tmp_dir, "adapters.fasta")
        write_adapters_fasta(args.adapter1, args.adapter2, adapters_fasta)
        
        vsearch_results = iterative_vsearch(
            input_bam=args.input_bam,
            adapters_fasta=adapters_fasta,
            output_dir=args.output_dir,
            tmp_dir=tmp_dir,
            min_len=args.min_len,
            id=args.identity,
            rounds=args.rounds,
            threads=args.threads
        )
        adapter_summary = vsearch_results.group_by("target").agg(pl.count().alias("Occurrences"))\
            .with_columns((pl.col("Occurrences") / (pl.col("Occurrences").sum()))\
                .alias("Frequency")).sort("Frequency", descending=True)
        print(adapter_summary)
    

    valid_pairs_dict, stat = find_valid_pairs(vsearch_results)
    config_stat = os.path.join(args.output_dir, "config_stat.tsv")
    with open(config_stat, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for key, value in stat.items():
            writer.writerow([key, value])
    
    outbam = os.path.join(args.output_dir,"extracted_reads.bam")
    concurrent_extract_reads(
        input_bam=args.input_bam,
        valid_pairs_dict=valid_pairs_dict,
        output_bam=outbam,
        threads=args.threads,
        batch_size = args.batch_size
    )

    print("Process completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads from basecalled bam file and collect statistics.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file path")
    parser.add_argument("-a1", "--adapter1", required=True, help="Sequence of the first adapter on the 5' end")
    parser.add_argument("-a2", "--adapter2", required=True, help="Sequence of the second adapter on the 3' end")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-ml", "--min_len", type=int, default=20, help="Minimum length of reads to be considered by vsearch")
    parser.add_argument("-id", "--identity", type=float, default=0.7, help="Minimum identity for adapter scanning")
    parser.add_argument("-r", "--rounds", type=int, default=3, help="Number of iterative rounds for vsearch")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("-s", "--batch_size", type=int, default=1000, help="Batch size for writing output files")
    args = parser.parse_args()
    main(args)

