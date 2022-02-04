import pandas as pd
from pandarallel import pandarallel
import argparse

# FUNCTIONS ===================================================================
# start with ATG ends with TAA, TGA, or TAG
def findORFSeq(seq, start="ATG", stop=["TAA", "TGA", "TAG"]):
    """
    find sequences starting with ATG and ending with termination codons
    """
    # look for atg
    s = None
    for i in enumerate(0, len(seq)):
        _3mer_s = seq[i : i + 3]
        if _3mer_s == start:
            s = i
            break
    f = None
    for j in enumerate(s, len(seq)):
        _3mer_f = seq[j : j + 3]
        if _3mer_f in stop:
            f = j + 3
    if (s is None) or (f is None):
        return None
    else:
        return seq[s:f]


def findString(seq, to_find="ATG"):
    """
    find a given substring on a longer string
    """
    # sanity check
    assert len(to_find) < len(seq), "ERROR: substring longer than sequence"

    # look for string on the sequence
    found_idx_lst = []
    len_str = len(to_find)
    seq_size = len(seq)
    for i in range(0, seq_size):
        if (i + len_str) > seq_size:
            break
        _mer_s = seq[i : i + len_str]
        if _mer_s == to_find:
            found_idx_lst.append(i)
    # return indexes of substring
    return found_idx_lst


def parseRandomFasta(in_file):
    """
    parse random fasta file.
    PS: it assumes sequence is on a single line

    return
    ------
    <dct>
    """
    fasta_in = open(in_file, "r")
    dct_lst = []

    for line in fasta_in:
        if line.startswith(">"):
            dct = {">": int(line[1 : len(line)].replace("\n", "").replace(" ", ""))}
            continue
        else:
            dct["seq"] = line.replace("\n", "")
            dct_lst.append(dct)
            continue
    return dct_lst


# for pandas
def countIdxs(row, colnm):
    return len(row[colnm])


# -----------------------------------------------------------------------------
# === INPUT HANDLING ==========================================================
dsc = "count and identify a given subsequence on sequences of a fasta file"
parser = argparse.ArgumentParser(description=dsc)
# --- general -----------------------------------------------------------------
parser.add_argument("fastaFile", type=str, default=None, help="input fasta file")

parser.add_argument(
    "subString",
    type=str,
    help="substring to look for on sequences (ex. 'ATG')",
    default=None,
)
parser.add_argument("-ncpus", type=int, default=1, help="number of cpus to use")

args = parser.parse_args()
# -----------------------------------------------------------------------------
# === MAIN ====================================================================
print("|--> INPUT")
print(f"| fastaFile = {args.fastaFile}")
print(f"| subString = {args.subString}")
print(f"| ncpus = {args.ncpus}")

# start pandas parallelilzation
pandarallel.initialize(nb_workers=args.ncpus)
# mount sequence dataframe
print("@ parsing fasta file")
seqs_dct = parseRandomFasta(args.fastaFile)
seqs_df = pd.DataFrame(seqs_dct)
new_col_nm = args.subString + "_pos"
print("@ searching for strings")
seqs_df[new_col_nm] = seqs_df["seq"].parallel_apply(findString, to_find=args.subString)
new_count_col = args.subString + "_count"
seqs_df[new_count_col] = seqs_df.parallel_apply(countIdxs, colnm=new_col_nm, axis=1)
print("|--| summary |--|")
print(f"| TOTAL FOUND = {seqs_df[new_count_col].sum()}")
print(f"| AVERAGE = {seqs_df[new_count_col].mean()} / seqs")
print("> saving results as '.csv' at current working dir...")
seqs_df[[">", new_col_nm, new_count_col]].to_csv("./results.csv")
print(":: DONE ::")
# TODO add find gene function
#  --> given a ATG index, loook for stop codons ahead of index, check if is
#    divisible by three, store sequence
