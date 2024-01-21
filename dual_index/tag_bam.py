import argparse
import pysam 
import pandas as pd

def add_tags(tags_file, in_bam, out_bam):
    """Add all the required tags to the BAM file."""
    
    df = pd.read_csv(tags_file, sep='\t', index_col="read_id")
    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in.fetch():
                read_id = align.query_name
                try:
                    row = df.loc[read_id]
                except KeyError:
                    continue  # No barcode/umi for this read
                # assigned well id
                align.set_tag('WI', row['wellid'], value_type="Z")
                # Corrected UMI = UB:Z
                align.set_tag("UB", row['UB'], value_type="Z")
                # Annotated gene name = GN:Z
                align.set_tag("GN", row['gene'], value_type="Z")


def main(args):
    add_tags(args.tags, args.in_bam, args.out_bam, args.chrom)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add tags to BAM file")
    parser.add_argument("-t", "--tags", required=True, help="Tag file")
    parser.add_argument("-i", "--in_bam", required=True, help="Input BAM file to be tagged")
    parser.add_argument("-o", "--out_bam", required=True, help="Output BAM file with tags")
    args = parser.parse_args()
    main(args)