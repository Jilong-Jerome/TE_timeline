import sys
sp_dict = {
        "dum":["lin","mim","bic","sar","ten"],
        "ten":["lin","mim","bic","sar","dum"],
        "sar":["lin","mim","ten","dum","bic"],
        "bic":["lin","mim","ten","dum","sar"]
        }
import random
from Bio import SeqIO

def find_closest_species_and_select_sequence(fasta_file, focus_species, sp_dict):
    species_in_file = set()
    sequences_by_species = {}

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        parts = record.id.split("|")
        species_id = parts[1].split("_")[0]
        species_in_file.add(species_id)
        if species_id not in sequences_by_species:
            sequences_by_species[species_id] = []
        sequences_by_species[species_id].append(record.id)

    # Find the closest species
    closest_species = focus_species
    for species in sp_dict[focus_species]:
        if species in species_in_file:
            closest_species = species
            break

    # Randomly select a sequence from the closest species
    selected_sequence = random.choice(sequences_by_species[closest_species])

    # Determine root situation
    root_situation = "hard to root" if closest_species in sp_dict[focus_species][-1:] else "rooted with outgroup"

    # Output data
    node_id = selected_sequence.split("|")[0]
    return [node_id, focus_species, closest_species, selected_sequence, root_situation]

# Example usage
fasta_file = sys.argv[1]  # Replace with your actual file path
focus_species = sys.argv[2]# Replace with your actual focus species
outname = sys.argv[3]
output_data = find_closest_species_and_select_sequence(fasta_file, focus_species, sp_dict)

# Writing the result to a file without using csv.writer
output_file = "{outname}_root.tsv".format(outname = outname)
with open(output_file, 'w') as file:
    # Joining each row's elements with a tab
    header = '\t'.join(['Node ID', 'Focus Species', 'Closest Species', 'Selected Sequence', 'Root Situation']) + '\n'
    file.write(header)
    data_line = '\t'.join(output_data) + '\n'
    file.write(data_line)

