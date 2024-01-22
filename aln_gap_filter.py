import sys

def read_fasta(filename):
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(filename, 'r') as file:
        sequence_id = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                sequence_id = line[1:]
                sequences[sequence_id] = ''
            else:
                sequences[sequence_id] += line
    return sequences

def filter_alignment(sequences, gap_threshold=0.95):
    """Filter alignment based on gap fraction at each position."""
    num_sequences = len(sequences)
    alignment_length = len(next(iter(sequences.values())))
    positions_to_keep = []

    for i in range(alignment_length):
        gap_count = sum(seq[i] == '-' for seq in sequences.values())
        if gap_count / num_sequences <= gap_threshold:
            positions_to_keep.append(i)

    filtered_sequences = {seq_id: ''.join(seq[i] for i in positions_to_keep)
                          for seq_id, seq in sequences.items()}
    return filtered_sequences

def write_fasta(sequences, filename):
    """Write sequences to a FASTA file."""
    with open(filename, 'w') as file:
        for seq_id, seq in sequences.items():
            file.write(f">{seq_id}\n")
            file.write(f"{seq}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]
    sequences = read_fasta(input_file)
    filtered_sequences = filter_alignment(sequences)
    write_fasta(filtered_sequences, output_file)

if __name__ == "__main__":
    main()

