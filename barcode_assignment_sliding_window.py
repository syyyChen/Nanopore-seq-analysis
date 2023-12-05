import os
import sys
import argparse
import csv
from Bio.Seq import Seq
import rapidfuzz
from rapidfuzz.process import extractOne


def get_whitelist(whitelist_file):
    """Read the whitelist file as dictionaries
    :param whitelist_file: path to the whitelist file in tsv format
    :return: two dictionaries containing the labels and indexes (one for reverse complementary strands)
    """
    whitelist_dict = {}
    whitelist_rc_dict = {}
    with open(whitelist_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            label, index = parts[0], parts[1]
            whitelist_dict[label] = index
            index_rc = str(Seq(index).reverse_complement())
            whitelist_rc_dict[label] = index_rc
    return whitelist_dict, whitelist_rc_dict



def barcode_scan(sequence, whitelist_dict, max_ed=2, min_ed_diff=1, flank=5, bc_length=10):
    """
    Scans a sequence against a whitelist of barcodes to find the best matching barcode.
    Args:
    sequence (str): The DNA sequence to be scanned.
    whitelist_dict (dict): A dictionary of barcodes with their indices.
    max_ed (int): Maximum allowed edit distance.
    min_ed_diff (int): Minimum edit distance difference to distinguish the top two candidates.
    flank (int): Number of flanking nucleotides on each side of the barcode.

    Returns:
    tuple: A tuple containing the best matching barcode, its edit distance, index, and window shift.
           Returns None for all elements if no suitable match is found.
    """
    sld_windows = [sequence[i:i+bc_length] for i in range(len(sequence) - bc_length + 1)] if len(sequence) > 10 else [sequence]
    candidates = []
    for idx, bc in whitelist_dict.items():
        res = extractOne(
            bc,
            sld_windows,
            scorer=rapidfuzz.distance.Levenshtein.distance
        )
        # keep the barcode idx, min_edit_distance, window_shift
        candidates.append((idx, res[1], res[2]-flank))
    
    candidates.sort(key=lambda x: x[1])
    top_hit_bc = whitelist_dict[candidates[0][0]]
    edit_distance = candidates[0][1]
    frameshift = candidates[0][2]
    
    if edit_distance <= max_ed:
        if candidates[1][1] - edit_distance >= min_ed_diff:
            bc_idx = candidates[0][0]
        else:
            bc_idx = 'Ambiguous'
    else:
        bc_idx = 'NoHit'

    return (top_hit_bc, edit_distance, bc_idx, frameshift)



def process_barcodes(uncorr_bc_file, output_dir, i7_whitelist, i5_whitelist, 
                     chunk_size=5000, max_ed=2, min_ed_diff=1, flank=5, bc_length=10):
    """
    Process a TSV file in chunks.
    Args:
    uncorr_bc_file(str): Path to the input TSV file with uncorrected barcodes.
    output_dir(str): Path to the output directory.
    i7_whitelist(str): Path to the i7 whitelist file
    i5_whitelist(str): Path to the i5 whitelist file
    chunk_size(int): Number of lines to read in each chunk.
    """
    corr_bc_file = os.path.join(output_dir, 'correct_barcodes.tsv')
    bc_counts_file = os.path.join(output_dir, 'readcounts_by_barcodes.tsv')
    
    (i7, i7rc) = get_whitelist(whitelist_file = i7_whitelist)
    (i5, i5rc) = get_whitelist(whitelist_file = i5_whitelist)
    wl_p = {'i7': i7, 'i5': i5rc}
    wl_n = {'i7': i7rc, 'i5': i5}
    bc_counts = {}

    
    with open(uncorr_bc_file, 'r') as input_file, open(corr_bc_file, 'w', newline='') as output_file:
        batch = []
        reader = csv.reader(input_file, delimiter='\t')
        column_names = next(reader)
        writer = csv.DictWriter(output_file, 
                                fieldnames=["read_id", "i7index_uncorr", "i7index_qual", "i7index_idx","i7index_ed", "i7index_corr", "i7index_frameshift", 
                                            "i5index_uncorr", "i5index_qual", "i5index_idx", "i5index_ed", "i5index_corr", "i5index_frameshift", "wellid"], 
                                delimiter='\t')
        writer.writeheader()


        for bc_record in reader:
            id, pattern, i7_uncorr, i7_qual, i5_uncorr, i5_qual = bc_record
            wl_dic = wl_p if pattern == '+' else wl_n
            (i7_bc, i7_ed, i7_idx, i7_shift) = barcode_scan(i7_uncorr, wl_dic['i7'], max_ed, min_ed_diff, flank, bc_length)
            (i5_bc, i5_ed, i5_idx, i5_shift) = barcode_scan(i5_uncorr, wl_dic['i5'], max_ed, min_ed_diff, flank, bc_length)
            wellid = f"{i5_idx}{i7_idx}" if i7_idx not in ['Ambiguous', 'NoHit'] and i5_idx not in ['Ambiguous', 'NoHit'] else 'Unsure'
            
            bc_corr = {
                'read_id': id,
                'i7index_uncorr': i7_uncorr,
                'i7index_qual': i7_qual,
                'i7index_idx': i7_idx,
                'i7index_ed': i7_ed,
                'i7index_corr': i7_bc,
                'i7index_frameshift': i7_shift,
                'i5index_uncorr': i5_uncorr,
                'i5index_qual': i5_qual,
                'i5index_idx': i5_idx,
                'i5index_ed': i5_ed,
                'i5index_corr': i5_bc,
                'i5index_frameshift': i5_shift,
                'wellid': wellid
            }                  
            
            bc_counts[wellid] = bc_counts.get(wellid, 0) + 1
            batch.append(bc_corr)

            if len(batch) == chunk_size:
                writer.writerows(batch)
                batch.clear()

        if batch:
            writer.writerows(batch)
        with open(bc_counts_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['wellid', 'read_count'])
            for wellid, count in bc_counts.items():
                writer.writerow([wellid, count])



def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    process_barcodes(
        uncorr_bc_file=args.input_bc_tsv,
        output_dir=args.output_dir,
        i7_whitelist=args.i7index,
        i5_whitelist=args.i5index,
        max_ed=args.max_ed,
        min_ed_diff=args.min_ed_diff,
        flank=args.flank,
        bc_length=args.bc_length
    )
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barcodes correction and assignment")
    parser.add_argument("-i", "--input_bc_tsv", required=True, help="Input tsv file recording barcodes")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-i7", "--i7index", required=True, help="Whitelist of i7 indexes")
    parser.add_argument("-i5", "--i5index", required=True, help="Whitelist of i5 indexes")
    parser.add_argument("-ed", "--max_ed", type=int, default=2, help="Maximum edit distance for barcode correction")
    parser.add_argument("-di", "--min_ed_diff", type=int, default=1, help="Minimum edit distance between first hits and second hits")
    parser.add_argument("-f", "--flank", type=int, default=5, help="Extracted extra bases flanking the barcode position")
    parser.add_argument("-l", "--bc_length", type=int, default=10, help="barcode length")
    args = parser.parse_args()
    main(args)