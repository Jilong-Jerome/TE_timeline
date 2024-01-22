import sys
def determine_condition(species_found, target_species, sister_species):
    target_present = target_species in species_found
    sister_present = sister_species in species_found

    if target_present and not sister_present and len(species_found) == 1:
        return 'A'  # Only target species present
    elif target_present and sister_present and len(species_found) == 2:
        return 'B'  # Both target and sister species present, no others
    elif len(species_found) >= 3:
        return 'C'  # Three or more species, including others
    else:
        return 'D'  # Any other situation

def process_node_data(node_data, target_species, sister_species, species_of_interest):
    species_count = {species: 0 for species in species_of_interest}
    species_found = set()

    for chromosome_id in node_data:
        species_code = chromosome_id.split('_')[0]
        if species_code in species_count:
            species_count[species_code] += 1
            species_found.add(species_code)

    condition = determine_condition(species_found, target_species, sister_species)
    return condition, species_count

def count_species_codes(input_file, target_species, sister_species, species_of_interest):
    node_data = {}
    
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            node_id, chromosome_id = parts[0], parts[1]

            if node_id not in node_data:
                node_data[node_id] = []
            node_data[node_id].append(chromosome_id)

    output_lines = []

    for node_id, data in node_data.items():
        condition, result = process_node_data(data, target_species, sister_species, species_of_interest)
        line = f"{node_id}\t{condition}"
        for species in species_of_interest:
            line += f"\t{result[species]}"
        output_lines.append(line)

    return output_lines



# Example usage
target_species_dict = {'dum': 'ten',
                       'ten': 'dum',
                       'sar':'bic',
                       'bic':'sar'}
species_of_interest = ['dum', 'ten', 'sar','bic' , 'mim','lin']  # Replace with actual species codes
input_file = sys.argv[1]
output_name = sys.argv[2]
target_species = sys.argv[3]
sister_species = target_species_dict[target_species]
output_lines = count_species_codes(input_file, target_species, sister_species, species_of_interest)

for line in output_lines:
    print(line)

# Optionally, write these lines to a file
with open(output_name, "w") as out_file:
    for line in output_lines:
        out_file.write(line + "\n")

