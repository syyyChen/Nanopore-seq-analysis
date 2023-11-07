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

valid_patterns = {
    "adapter1_f": "adapter2_f",
    "adapter2_r": "adapter1_r"
}

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


def mask_regions(vsearch_res_df, fasta, tmp_dir):
    
    regions_to_mask = defaultdict(list)
    for row in vsearch_res_df.iter_rows():
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


def iterative_vsearch(input_bam, adapter1, adapter2, output_dir, min_len=20, id=0.7, rounds=3, threads=1):
    output = os.path.join(output_dir, "adapter_scan.vsearch.tsv")
    res_all = pl.DataFrame(schema=schema)
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        adapters_fasta = os.path.join(tmp_dir, "adapters.fasta")
        write_adapters_fasta(adapter1, adapter2, adapters_fasta)
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
    valid_pairs = {}
    for query_id, grouped_df in vs_result.group_by('query'):
        grouped_df = grouped_df.sort('qilo')
        if grouped_df['target'][0] == '*':
            continue
        for i in range(grouped_df.height - 1):
            current_adapter = grouped_df.row(i)
            next_adapter = grouped_df.row(i + 1)
            if (current_adapter[1] in valid_patterns and
                valid_patterns[current_adapter[1]] == next_adapter[1] and
                current_adapter[7] < next_adapter[6]):
                valid_pairs.setdefault(query_id, [])
                valid_pairs[query_id].append((current_adapter[7], next_adapter[6] - 1))
    return valid_pairs


def extract_reads(input_bam, valid_pairs_dict, output_dir):
    
    output_bam = os.path.join(output_dir, "extracted_reads.bam")
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bamfile:
        
        with pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as outfile:
            for read in bamfile.fetch(until_eof=True):
                query_name = read.query_name
                if query_name in valid_pairs_dict:
                    for idx, (start, end) in enumerate(valid_pairs_dict[query_name]):
                        new_read = pysam.AlignedSegment()
                        new_read.query_name = f"{query_name}_{idx}"
                        new_read.flag = read.flag
                        new_read.seq = read.seq[start:end]  
                        new_read.qual = read.qual[start:end] if read.qual else None
                        new_read.tags = read.tags 
                        outfile.write(new_read)


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
    plt.xticks(np.arange(min(lengths), max(lengths) + 200, 200))
    plt.axvline(median_length, color='red', linestyle='dashed', linewidth=1)
    plt.axvline(average_length, color='orange', linestyle='dashed', linewidth=1)
    plt.text(median_length, max(frequencies)*0.8, f'Median: {median_length:.2f}', color='red')
    plt.text(average_length, max(frequencies)*0.7, f'Average: {average_length:.2f}', color='orange')
    plt.grid(True)  

    output_path = os.path.join(output_dir, 'length_distribution.png')
    plt.savefig(output_path)
    plt.close()  


def main(args):
    
    vsearch_results = iterative_vsearch(
        input_bam=args.input_bam,
        adapter1=args.adapter1,
        adapter2=args.adapter2,
        output_dir=args.output_dir,
        min_len=args.min_len,
        id=args.identity,
        rounds=args.rounds,
        threads=args.threads
    )

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
    parser.add_argument("-a1", "--adapter1", required=True, help="Sequence of the first adapter")
    parser.add_argument("-a2", "--adapter2", required=True, help="Sequence of the second adapter")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-ml", "--min_len", type=int, default=20, help="Minimum length of reads to consider by vsearch")
    parser.add_argument("-id", "--identity", type=float, default=0.7, help="Minimum identity for adapter scanning")
    parser.add_argument("-r", "--rounds", type=int, default=3, help="Number of iterative rounds for vsearch")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for vsearch")

    args = parser.parse_args()
    main(args)





