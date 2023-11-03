import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from collections import defaultdict

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

def read_vsearch_res(tsv_file, selected_columns=vsearch_colnames):
    
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
    

def iterative_vsearch(input_fasta, adapters_fasta, output, tmp_dir, id=0.7, rounds=3, threads=1):
    
    query = input_fasta
    tmp_fasta = os.path.join(tmp_dir, "tmp.fasta")
    res_all = pl.DataFrame(schema=schema)
    
    for round in range(1, rounds+1):
        result = os.path.join(tmp_dir, "vsearch_res_{no}.tsv".format(no=round))
        vsearch_cmd = "vsearch --usearch_global {query_fasta} --db {adapters} --minseqlength 20 --maxaccepts 5 --id {min_id} --strand plus \
        --wordlength 3 --minwordmatches 10 --threads {nthreads} --userfields \
        'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
        --userout {out} --matched {tmp_file}".format(
            query_fasta=query,
            adapters=adapters_fasta,
            min_id=id,
            nthreads=threads,
            out=result,
            tmp_file=tmp_fasta
        )
        if round == 1:
            vsearch_cmd += " --output_no_hits"
            
        p = subprocess.Popen(vsearch_cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        p.communicate()
        
        # stop the iteration if there are no more targets found
        if not os.path.getsize(result):
            break
        vsearch_results = read_vsearch_res(result)
        res_all = pl.concat([res_all, vsearch_results])
        # reads found at least one target will be used for the next round of searching
        # targets already found are masked with U
        if round < rounds:
            mask_regions(vsearch_results, tmp_fasta, tmp_dir)
            subprocess.run("rm {tmp_fasta}".format(tmp_fasta=tmp_fasta), shell=True)
            query = os.path.join(tmp_dir, "masked.fasta")
        
    subprocess.run("rm -r {tmp}/*".format(tmp=tmp_dir), shell=True)
    res_all=res_all.sort(["query", "qilo"])
    res_all.write_csv(output, separator="\t", has_header=False)
    