from Bio import SeqIO
from Bio.Seq import Seq
import sys

def load_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences[record.id] = record.seq
    return sequences

def get_strand_direction_and_sequence(subject_seq, start, end):
    if start <= end:
        return '+', subject_seq[start - 1:end]  # 1-based to 0-based indexing for start
    else:
        # Reverse complement for reverse strand
        return '-', subject_seq[end - 1:start].reverse_complement() 

def extract_subsequences(blast_file, fasta_sequences, output_file):
    with open(blast_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split()
            node = fields[0]
            subject_id = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            length = fields[4]

            strand, subsequence = get_strand_direction_and_sequence(fasta_sequences[subject_id], start, end)
            header = f">{node}|{subject_id}|{start}|{end}|{strand}|{length}"
            outfile.write(f"{header}\n{subsequence}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_subsequences.py <blast_output.txt> <genome_fasta.fa> <output_fasta.fa>")
        sys.exit(1)

    blast_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    fasta_sequences = load_fasta(fasta_file)
    extract_subsequences(blast_file, fasta_sequences, output_file)

if __name__ == "__main__":
    main()

