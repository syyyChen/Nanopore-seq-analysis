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


def find_valid_pairs(vs_result):
    print("Finding valid pairs...")
    stat = {}
    stat['no_adapter'] = vs_result.filter(pl.col('target') == '*').height
    n_with_adapter = vs_result.filter(pl.col('target') != '*').n_unique(subset = ['query'])
    
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
        i7 = row[4]
        i5 = row[5]
        valid_pairs_dict[query_id] = list(zip(start_positions, end_positions, patterns, i7, i5))
    
    print(f"Found {n_valid} reads with valid pairs.")
    return valid_pairs_dict, stat


def valid_stat_plot(data_dict, output_dir):
    def custom_autopct(pct, allvals):
        absolute = int(pct/100.*sum(allvals))
        return "{:.1f}%({:d})".format(pct, absolute)
    labels = list(data_dict.keys())
    sizes = list(data_dict.values())
    wedgeprops = {"edgecolor": "white", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True}
    plt.figure(figsize=(12, 6)) 
    fig, ax = plt.subplots()
    wedges, _, autotexts = ax.pie(sizes, autopct=lambda pct: custom_autopct(pct, sizes), wedgeprops=wedgeprops)
    for autotext in autotexts:
        autotext.set_size(8)
    plt.legend(wedges, labels, loc='upper left', bbox_to_anchor=(1.2, 1.1))
    output_path = os.path.join(output_dir, 'config_stat.png')
    plt.savefig(output_path)
    plt.close()  


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
                    multi_reads = 1 if len(valid_pairs_dict[query_name])> 1 else 0
                    for idx, (start, end, pattern, i7_start, i5_start) in enumerate(valid_pairs_dict[query_name]):
                        i7 = read.seq[i7_start:i7_start+10]
                        i5 = read.seq[i5_start:i5_start+10]
                        new_read = pysam.AlignedSegment()
                        new_read.query_name = f"{query_name}_{idx}"
                        new_read.flag = read.flag
                        new_read.seq = read.seq[start:end]  
                        new_read.qual = read.qual[start:end]
                        new_read.tags = read.tags 
                        new_read.set_tag('mp', multi_reads, 'i')
                        new_read.set_tag('pn', pattern, 'Z')
                        new_read.set_tag('i7', i7, 'Z')
                        new_read.set_tag('i5', i5, 'Z')
                        umi = read.seq[end-8:end] if pattern=='+' else read.seq[start:start+8]
                        new_read.set_tag('mi', umi, 'Z')
                        
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
        for start, end, pattern, i7, i5 in values:
            lengths.append(end - start)

    median_length = np.median(lengths)
    average_length = np.mean(lengths)

    counts, bin_edges = np.histogram(lengths, bins=50)
    total_counts = sum(counts)
    frequencies = counts / total_counts

    plt.figure(figsize=(15, 9))  
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
    

    valid_pairs_dict, stat = find_valid_pairs(vsearch_results)
    
    valid_stat_plot(stat, output_dir=args.output_dir)

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





