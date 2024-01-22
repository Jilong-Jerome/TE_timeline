import sys
from Bio import SeqIO

def read_sequence_names(filename):
    node = filename.split("/")[-1].split("_")[0]
    with open(filename, 'r') as file:
        return [node+"|"+line.strip() for line in file]

def filter_fasta(fasta_file, sequence_names, output_file):
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in sequence_names:
                SeqIO.write(record, outfile, 'fasta')

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <fasta_file> <sequence_names_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence_names_file = sys.argv[2]
    output_file = sys.argv[3]

    sequence_names = read_sequence_names(sequence_names_file)
    print(sequence_names)
    filter_fasta(fasta_file, sequence_names, output_file)

if __name__ == "__main__":
    main()

