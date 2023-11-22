import os
import sys
import argparse
import csv
from Bio.Seq import Seq
import rapidfuzz
from rapidfuzz.process import extract



# import the whitelist
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


# process single barcode and return results
def calc_ed_with_whitelist(bc_uncorr, whitelist_dict, score_cutoff=3, max_ed=2, min_ed_diff=1):
    """Find barcodes in a whilelist with the smallest edit distance.
    :param bc_uncorr: uncorrected barcode
    :param whitelist_dict: dictionary of possible barcodes
    :param: score_cutoff: rapidfuzz score cutoff - edit distances above this are not reported.
    :return:
        best matching barcode
        label for the best matching barcode
        edit distance between best match and uncorrected barcode
        edit distance difference between top match and second top match
    """
    
    result = extract(
        bc_uncorr,
        whitelist_dict,
        scorer=rapidfuzz.distance.Levenshtein.distance,
        score_cutoff=score_cutoff)
    
    corr_condition = (
        ((len(result) == 1) and (result[0][1] <= max_ed)) or
        ((len(result) > 1) and (result[0][1] <= max_ed) \
            and (result[1][1] - result[0][1] >= min_ed_diff))
    )
    
    if corr_condition:
        bc_match = result[0][0]
        bc_match_ed = result[0][1]
        bc_match_idx = result[0][2]
        #next_match_diff = result[1][1] - bc_match_ed if len(result) > 1 \
        #    else score_cutoff - bc_match_ed
    else:
        bc_match = None
        bc_match_ed = None
        bc_match_idx = None
        #next_match_diff = None
    
    return bc_match, bc_match_ed, bc_match_idx


def process_barcodes(uncorr_bc_file, output_dir, i7_whitelist, i5_whitelist, chunk_size=5000, score_cutoff=3, max_ed=2, min_ed_diff=1):
    """
    Process a TSV file in chunks.
    :param uncorr_bc_file: Path to the input TSV file with uncorrected barcodes.
    :param output_dir: Path to the output directory.
    :param i7_whitelist: Path to the i7 whitelist file
    :param i5_whitelist: Path to the i5 whitelist file
    :param chunk_size: Number of lines to read in each chunk.
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
                                fieldnames=["read_id", "i7index_uncorr", "i5index_uncorr", \
                                    "i7index_ed", "i5index_ed", "i7index_corr", "i5index_corr", "wellid"], delimiter='\t')
        writer.writeheader()


        for bc_record in reader:
            id, pattern, i7_uncorr, i5_uncorr = bc_record
            wl_dic = wl_p if pattern == '+' else wl_n
            (i7_bc, i7_ed, i7_idx) = calc_ed_with_whitelist(i7_uncorr, wl_dic['i7'], score_cutoff, max_ed, min_ed_diff)
            (i5_bc, i5_ed, i5_idx) = calc_ed_with_whitelist(i5_uncorr, wl_dic['i5'], score_cutoff, max_ed, min_ed_diff)
            wellid = f"{i5_idx}{i7_idx}" if i7_idx is not None and i5_idx is not None else 'Unsure'
            
            bc_corr = {
                'read_id': id,
                'i7index_uncorr': i7_uncorr,
                'i5index_uncorr': i5_uncorr,
                'i7index_ed': i7_ed,
                'i5index_ed': i5_ed,
                'i7index_corr': i7_bc,
                'i5index_corr': i5_bc,
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
        score_cutoff=args.cutoff
    )
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barcodes correction and assignment")
    parser.add_argument("-i", "--input_bc_tsv", required=True, help="Input tsv file recording barcodes")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-i7", "--i7index", required=True, help="Whitelist of i7 indexes")
    parser.add_argument("-i5", "--i5index", required=True, help="Whitelist of i5 indexes")
    parser.add_argument("-ed", "--max_ed", type=int, default=2, help="Maximum edit distance for barcode correction")
    parser.add_argument("-di", "--min_ed_diff", type=int, default=1, help="Minimum edit distance between first hits and second hits")
    parser.add_argument("-c", "--cutoff", type=int, default=3, help="Cutoff edit distance for whitelist searching")
    args = parser.parse_args()
    main(args)

