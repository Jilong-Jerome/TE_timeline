import sys
from Bio import AlignIO
import numpy as np
import pandas as pd

def calculate_distance(seq1, seq2, gap_open_penalty, gap_extend_penalty):
    """Calculate distance as scaled mismatch count with gap penalties."""
    mismatch_count, gap_open, comparable_length = 0, False, 0

    for a, b in zip(seq1, seq2):
        if a != '-' or b != '-':
            comparable_length += 1
            if a == '-' or b == '-':
                if not gap_open:
                    mismatch_count += gap_open_penalty
                    gap_open = True
                else:
                    mismatch_count += gap_extend_penalty
            else:
                gap_open = False
                if a != b:
                    mismatch_count += 1

    if comparable_length > 0:
        scaled_mismatch = mismatch_count / comparable_length
    else:
        scaled_mismatch = -1  # Denote missing values with -1

    return scaled_mismatch

def main():
    if len(sys.argv) != 5:
        print("Usage: python script.py <alignment_file.fasta> <gap_open_penalty> <gap_extend_penalty> <output_distance>")
        sys.exit(1)

    alignment_file = sys.argv[1]
    gap_open_penalty = float(sys.argv[2])
    gap_extend_penalty = float(sys.argv[3])
    out_dist_name = sys.argv[4]
    # Read the MSA file
    alignment = AlignIO.read(alignment_file, "fasta")

    # Initialize the distance matrix
    num_seqs = len(alignment)
    distance_matrix = np.zeros((num_seqs, num_seqs))

    # Calculate distances
    for i in range(num_seqs):
        for j in range(i+1, num_seqs):
            distance = calculate_distance(alignment[i].seq, alignment[j].seq, gap_open_penalty, gap_extend_penalty)
            distance_matrix[i][j] = distance_matrix[j][i] = distance

    # Convert to DataFrame for easy CSV output
    seq_ids = [record.id for record in alignment]
    matrix_df = pd.DataFrame(distance_matrix, index=seq_ids, columns=seq_ids)
    matrix_df.to_csv(out_dist_name)

if __name__ == "__main__":
    main()

