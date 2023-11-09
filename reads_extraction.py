import os
import argparse
import sys
import numpy as np
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from collections import defaultdict
import tempfile
import pysam
import matplotlib.pyplot as plt

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


def mask_regions(vs_result_df, fasta, tmp_dir):
    """Mask mapped regions of reads with 'U' according to vsearch results """
    regions_to_mask = defaultdict(list)
    for row in vs_result_df.iter_rows():
        if row[11] == 0:
            continue
        query_id = row[0]
        start = row[6] - 1
        end = row[7]
        regions_to_mask[query_id].append((start, end))
    
    masked_records = []
    for record in SeqIO.parse(fasta, "fasta"):
        masked_seq_str = str(record.seq)
        for start, end in regions_to_mask[record.id]:
            masked_seq_str = masked_seq_str[:start] + 'U' * (end - start) + masked_seq_str[end:]
        masked_record = SeqRecord(Seq(masked_seq_str), id=record.id, description=record.description)
        masked_records.append(masked_record)
    SeqIO.write(masked_records, os.path.join(tmp_dir, "masked.fasta"), "fasta")


def iterative_vsearch(input_bam, adapters_fasta, output_dir, tmp_dir, min_len=20, id=0.7, rounds=3, threads=1):
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


def collect_config_statistics(vs_result, output_dir):
    """Collect statistics of adapter patterns in one read"""
    grouped_res = vs_result.select(['query','target']).group_by('query').agg(
        pl.col("target").str.concat("-").alias("Configuration"))

    # Calculate the occurrences and frequencies
    config_counts = grouped_res.group_by("Configuration").agg(
        pl.count().alias("Occurrences")
    ).with_columns(Frequency=(pl.col.Occurrences / grouped_res.height) * 100)

    config_counts = config_counts.sort("Frequency", descending=True)
    config_counts.write_csv(os.path.join(output_dir, 'config_stat.csv'))


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
        
    # Filter for valid pairs
    valid_pairs_condition = (
        ((vs_result['target'] == valid_patterns[0][0]) & (vs_result['next_target'] == valid_patterns[0][1])) |
        ((vs_result['target'] == valid_patterns[1][0]) & (vs_result['next_target'] == valid_patterns[1][1]))
    ) & (vs_result['qihi'] < vs_result['next_qilo'] - 1)
    valid_pairs_df = vs_result.filter(valid_pairs_condition)
    
    n_adapters = vs_result.filter(pl.col('target') != '*').height
    n_valid = valid_pairs_df.height
    n_unpaired = n_adapters - 2 * n_valid
    
    # Group and aggregate
    grouped_valid_pairs = valid_pairs_df.group_by('query').agg([
        pl.col('qihi').alias('start_positions'),
        pl.col('next_qilo').alias('end_positions')
    ])

    # Create dictionary of query to valid pairs
    valid_pairs_dict = {}
    for row in grouped_valid_pairs.iter_rows():
        query_id = row[0]
        start_positions = row[1]
        end_positions = [x - 1 for x in row[2]]
        valid_pairs_dict[query_id] = list(zip(start_positions, end_positions))
    
    print(f"Found {n_valid} valid pairs. {n_unpaired} adapters cannot be paired.")
    return valid_pairs_dict


def extract_reads(input_bam, valid_pairs_dict, output_dir, batch_size=1000):
    """Extract reads from the bam file with batch writing."""
    print("Extracting reads...")
    n_reads = 0
    output_bam = os.path.join(output_dir, "extracted_reads.bam")

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bamfile:    
        with pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as outfile:
            batch = []
            for read in bamfile.fetch(until_eof=True):
                query_name = read.query_name
                if query_name in valid_pairs_dict:
                    for idx, (start, end) in enumerate(valid_pairs_dict[query_name]):
                        new_read = pysam.AlignedSegment()
                        new_read.query_name = f"{query_name}_{idx}"
                        new_read.flag = read.flag
                        new_read.seq = read.seq[start:end]  
                        new_read.qual = read.qual[start:end]
                        new_read.tags = read.tags 
                        batch.append(new_read)
                        n_reads += 1

                        if len(batch) >= batch_size:
                            for b_read in batch:
                                outfile.write(b_read)
                            batch.clear()
            # Write remaining reads in the last batch
            for b_read in batch:
                outfile.write(b_read)

    print(f"{n_reads} reads with pattern have been extracted.")


def length_distribution_plot(valid_pairs_dict, output_dir):
    lengths = []
    for values in valid_pairs_dict.values():
        for start, end in values:
            lengths.append(end - start)

    median_length = np.median(lengths)
    average_length = np.mean(lengths)

    counts, bin_edges = np.histogram(lengths, bins=50)
    total_counts = sum(counts)
    frequencies = counts / total_counts

    plt.figure(figsize=(10, 6))  
    plt.bar(bin_edges[:-1], frequencies, width=np.diff(bin_edges), color='grey', alpha=0.7)  
    plt.title('Length Distribution of Extracted Reads')  
    plt.xlabel('Length')  
    plt.ylabel('Frequency') 
    stick_bin = max(lengths)//22//100*100
    plt.xticks(np.arange(0, max(lengths)+stick_bin, stick_bin))
    plt.axvline(median_length, color='red', linestyle='dashed', linewidth=1)
    plt.axvline(average_length, color='orange', linestyle='dashed', linewidth=1)
    plt.text(median_length, max(frequencies)*0.8, f'Median: {median_length:.2f}', color='red')
    plt.text(average_length, max(frequencies)*0.7, f'Average: {average_length:.2f}', color='orange')
    plt.grid(True)  

    output_path = os.path.join(output_dir, 'length_distribution.png')
    plt.savefig(output_path)
    plt.close()  


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
    
    collect_config_statistics(vsearch_results, args.output_dir)

    valid_pairs_dict = find_valid_pairs(vsearch_results)

    extract_reads(
        input_bam=args.input_bam,
        valid_pairs_dict=valid_pairs_dict,
        output_dir=args.output_dir
    )

    length_distribution_plot(
        valid_pairs_dict=valid_pairs_dict,
        output_dir=args.output_dir
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
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads of vsearch")
    args = parser.parse_args()
    main(args)





