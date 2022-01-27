import random
import argparse
import time

# -- FUNCTIONS ----------------------------------------------------------------
def random_string_generator(str_size, allowed_chars):
    return "".join(random.choice(allowed_chars) for x in range(str_size))


def generateRandomFasta(n_seqs, len_seqs, out_flnm):
    """
    generate a fasta file with random dna sequences
    """
    # generate file
    fasta_out = open(out_flnm, "w")
    # generate sequences and write it on the new file
    for i in range(0, n_seqs):
        # write header
        fasta_out.write(f"> {i} \n")
        # write seq
        random_seq = random_string_generator(len_seqs, ["A", "T", "C", "G"])
        fasta_out.write(random_seq + "\n")
    fasta_out.close()


# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # ---- INPUT --------------------------------------------------------------
    dsc = "generate a fasta file with a given number of sequences of a given size"
    parser = argparse.ArgumentParser(description=dsc)
    # --- general -------------------------------------------------------------
    parser.add_argument(
        "n_seqs", type=int, default=100, help="number of sequences to generate"
    )
    parser.add_argument(
        "len_seqs", type=int, default=200, help="size of sequences to generate"
    )

    parser.add_argument(
        "out_path", type=str, default="random.fasta", help="fasta output file name"
    )
    args = parser.parse_args()

    # generate fasta file
    print("@ generating random sequences...")
    start = time.time()
    generateRandomFasta(args.n_seqs, args.len_seqs, args.out_path)
    end = time.time()
    print(": total time = ", end - start)
    print(":: DONE ::")
