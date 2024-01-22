import sys
from Bio import SeqIO
from Bio.Seq import Seq

def translate_sequence(seq, frame):
    """Translate a DNA sequence from a given frame into an amino acid sequence."""
    return str(Seq(seq[frame:]).translate(table=1, to_stop=False, cds=False))

def find_orfs(sequence, min_len, max_len):
    """Find ORFs in a given sequence."""
    orfs = []
    for strand, nuc in [(+1, sequence), (-1, str(Seq(sequence).reverse_complement()))]:
        for frame in range(3):
            trans = translate_sequence(nuc, frame)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                if trans[aa_start] == 'M':
                    aa_end = trans.find('*', aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if (aa_end - aa_start) >= min_len and (aa_end - aa_start) <= max_len:
                        orf_seq = trans[aa_start:aa_end]
                        if strand == 1:
                            start = frame + aa_start * 3
                            end = frame + aa_end * 3
                        else:
                            start = len(sequence) - (frame + aa_end * 3)
                            end = len(sequence) - (frame + aa_start * 3)
                        orfs.append((start, end, strand, orf_seq))
                aa_start += 1
    return orfs

def format_orfs_to_fasta(orfs, prefix, original_seq_name):
    """Format ORFs to FASTA."""
    fasta_format = ""
    for i, (start, end, strand, orf_seq) in enumerate(sorted(orfs, key=lambda x: len(x[3]), reverse=True), 1):
        direction = '+' if strand == 1 else '-'
        header = f">{prefix}|ORF_{i}|{start + 1}|{end}|{direction}|{len(orf_seq)}|{original_seq_name}\n"
        fasta_format += header + orf_seq + "\n"
    return fasta_format

def get_prefix(header, field_number):
    """Extract the prefix from the sequence header based on the specified field number."""
    fields = header.split('|')
    try:
        return fields[field_number]
    except IndexError:
        print(f"Warning: Field number {field_number} is out of range for the header '{header}'. Using the full header as prefix.")
        return header

def display_help():
    """Display help message and usage instructions."""
    help_message = """
    ORF Finder Script
    -----------------
    This script finds all open reading frames (ORFs) in a DNA sequence provided in a FASTA file and outputs the results in a new FASTA file.

    Usage:
    python script.py <input_fasta_file> <output_fasta_file> <field_number>

    Arguments:
    <input_fasta_file> - Path to the input FASTA file containing DNA sequences.
    <output_fasta_file> - Path for the output FASTA file to write the ORFs.
    <field_number> - Field number in the FASTA header from which to extract the prefix for naming the ORFs.

    The output FASTA file will contain the identified ORFs in the following format:
    >{prefix}|ORF_{n}|{start_position}|{end_position}|{strand_direction}|{amino_acid_length}|{original_sequence_name}
    Each ORF sequence will be presented below its corresponding header.
    """
    print(help_message)

def main():
    if len(sys.argv) != 4:
        display_help()
        sys.exit(1)

    fasta_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    field_number = int(sys.argv[3])

    # Ensure the output file is empty before appending new content
    open(output_file_path, 'w').close()

    with open(fasta_file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            original_seq_name = record.id
            prefix = get_prefix(original_seq_name, field_number)
            orfs = find_orfs(sequence, min_len=30, max_len=float('inf'))  # ORF length in amino acids
            fasta_orfs = format_orfs_to_fasta(orfs, prefix, original_seq_name)

            with open(output_file_path, 'a') as output_file:
                output_file.write(fasta_orfs)

    print(f"ORFs written to {output_file_path}")

if __name__ == "__main__":
    main()

