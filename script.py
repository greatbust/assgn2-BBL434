# Nucleotide Sequence Alignment Tool
# Global Alignment
# Affine Gap Penalty (gap opening + gap extension)


import numpy as np

# fasta reader
def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()


# Alignment function
def Alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):

    n = len(seq1)
    m = len(seq2)

    Neg_inf = -10**9

    # Three matrices:
    # M = match/mismatch
    # X = gap in seq1
    # Y = gap in seq2

    M = np.full((n+1, m+1), Neg_inf)
    X = np.full((n+1, m+1), Neg_inf)
    Y = np.full((n+1, m+1), Neg_inf)

    traceback = {}

    M[0][0] = 0

    # Initialization
    for i in range(1, n+1):
        X[i][0] = gap_open + (i-1)*gap_extend

    for j in range(1, m+1):
        Y[0][j] = gap_open + (j-1)*gap_extend

    # Fill matrices
    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            M[i][j] = max(
                M[i-1][j-1],
                X[i-1][j-1],
                Y[i-1][j-1]
            ) + score

            X[i][j] = max(
                M[i-1][j] + gap_open,
                X[i-1][j] + gap_extend
            )

            Y[i][j] = max(
                M[i][j-1] + gap_open,
                Y[i][j-1] + gap_extend
            )

    # Final score
    i, j = n, m
    final_scores = [M[i][j], X[i][j], Y[i][j]]
    state = np.argmax(final_scores)

    # Traceback
    align1 = []
    align2 = []

    while i > 0 or j > 0:

        if state == 0:  # M matrix
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            if M[i][j] == M[i-1][j-1] + score:
                state = 0
            elif M[i][j] == X[i-1][j-1] + score:
                state = 1
            else:
                state = 2

            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1

        elif state == 1:  # X matrix (gap in seq2)
            if X[i][j] == M[i-1][j] + gap_open:
                state = 0
            else:
                state = 1

            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1

        elif state == 2:  # Y matrix (gap in seq1)
            if Y[i][j] == M[i][j-1] + gap_open:
                state = 0
            else:
                state = 2

            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    align1.reverse()
    align2.reverse()

    return "".join(align1), "".join(align2), max(final_scores)



if __name__ == "__main__":

    file1 = input("Enter FASTA file 1: ")
    file2 = input("Enter FASTA file 2: ")

    match = int(input("Match score: "))
    mismatch = int(input("Mismatch penalty: "))
    gap_open = int(input("Gap opening penalty: "))
    gap_extend = int(input("Gap extension penalty: "))

    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)

    aln1, aln2, score = Alignment(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    print("\n Best Alignment")
    print(aln1)
    print(aln2)
    print("\n Alignment Score:", score)
