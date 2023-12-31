import regex
import csv
import argparse
import sys


complement_trans = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn")


def find_best_match(pattern, sequence, max_ed, expected_start, expected_end):
    matches = regex.finditer(f"({pattern}){{e<={max_ed}}}", sequence, overlapped=True, flags=regex.ENHANCEMATCH)
    best_match = None
    best_score = 10

    for match in matches:
        start, end = match.start(), match.end()
        edit_distance = sum(match.fuzzy_counts)
        shift = max(abs(expected_start - start), abs(expected_end - end))
        score = edit_distance + shift
        if score < best_score:
            best_score = score
            best_match = match

    return best_match


def extract_umi(input_file, output_file, umi_len, ts_pattern, max_ed, flank, batch_size=5000):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, fieldnames=['read_id', 'pattern', 'umi', 'umi_qual'], delimiter='\t')
        writer.writeheader()

        batch = []
        for row in reader:
            strand = row['pattern']
            sequence = row['umi_uncorr']
            qual = row['umi_uncorr_qual']

            if strand == '+':
                sequence = sequence[::-1].translate(complement_trans)
                qual = qual[::-1]

            pattern = ts_pattern.replace('W', '[AT]')
            expected_start = flank + umi_len
            expected_end = flank + umi_len + len(ts_pattern)
            match = find_best_match(pattern, sequence, max_ed, expected_start, expected_end)

            if match:
                start, end = match.start(), match.end()
                umi = sequence[start - umi_len:start] if start>=umi_len else sequence[0:start]
                umi_qual = qual[start - umi_len:start] if start>=umi_len else qual[0:start]
            else:
                umi, umi_qual = '', ''

            batch.append({'read_id': row['read_id'], 'umi': umi, 'umi_qual': umi_qual})
            if len(batch) >= batch_size:
                writer.writerows(batch)
                batch.clear()
        
        if batch:
            writer.writerows(batch)


def main(args):
    extract_umi(
        input_file=args.input, 
        output_file=args.output, 
        umi_len=args.umi_len,
        ts_pattern=args.ts_pattern,
        max_ed=args.max_ed, 
        flank=args.flank
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequencing reads to extract UMIs.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with sequencing reads")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for extracted and orientated UMIs")
    parser.add_argument("-p", "--ts_pattern", type=str, required=True, help="Template switch oligo pattern")
    parser.add_argument("-ul", "--umi_len", type=int, default=8, help="Length of the UMI")
    parser.add_argument("-d", "--max_ed", type=int, default=1, help="Maximum edit distance allowed in the pattern oligos")
    parser.add_argument("-f", "--flank", type=int, default=5, help="Extra range flanking the expected pattern position")
    args = parser.parse_args()
    main(args)

