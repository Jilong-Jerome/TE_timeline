import sys
from Bio import SeqIO

def read_tab_file(tab_file):
    mapping = {}
    with open(tab_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            node, name = line.strip().split('\t')
            mapping[name] = node
    return mapping

def filter_and_rename_sequences(fasta_file, mapping, output_file):
    with open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract the name part after the second '|'
            name_part = record.id
            if name_part in mapping:
                new_name = f"{mapping[name_part]}|{name_part}"
                record.id = new_name
                record.description = new_name  # Update description as well
                SeqIO.write(record, out_fasta, "fasta")

def main():
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <tab_file> <fasta_file> <output_file>")
        print("\nThis script filters and renames sequences from a FASTA file based on a tab-separated file.")
        print("  <tab_file>: A tab-separated file with two columns: 'node' and 'name'.")
        print("  <fasta_file>: A FASTA file containing sequences to be filtered and renamed.")
        print("  <output_file>: The output FASTA file with filtered and renamed sequences.")
        return

    tab_file, fasta_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

    # Process the files
    mapping = read_tab_file(tab_file)
    filter_and_rename_sequences(fasta_file, mapping, output_file)

if __name__ == "__main__":
    main()

