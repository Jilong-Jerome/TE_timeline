import sys

def read_fasta(filename):
    """Read a FASTA file and return a dictionary of sequences."""
    with open(filename, 'r') as file:
        sequences = {}
        seq_name = ''
        seq_data = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_name:
                    sequences[seq_name] = seq_data
                seq_name = line[1:]  # Exclude the '>'
                seq_data = ''
            else:
                seq_data += line
        if seq_name:
            sequences[seq_name] = seq_data
    return sequences

def write_fasta(sequences, output_filename):
    """Write sequences to a FASTA file."""
    with open(output_filename, 'w') as file:
        for seq_name, seq_data in sequences.items():
            file.write(f">{seq_name}\n{seq_data}\n")

def deduplicate_sequences(input_filename, output_filename):
    """Read a FASTA file, deduplicate sequences, and write to a new file."""
    sequences = read_fasta(input_filename)
    write_fasta(sequences, output_filename)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    deduplicate_sequences(input_filename, output_filename)

if __name__ == '__main__':
    main()

