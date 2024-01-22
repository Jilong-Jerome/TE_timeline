import pickle
def process_genomic_file(file_path):
    genomic_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            node_id, chrom_id, start, end, length = line.strip().split('\t')
            if node_id not in genomic_data:
                genomic_data[node_id] = []
            genomic_data[node_id].append((chrom_id, int(start), int(end), int(length)))
    return genomic_data

def merge_genomic_data(file_paths):
    all_genomic_data = {}
    for file_path in file_paths:
        file_data = process_genomic_file(file_path)
        for node_id, regions in file_data.items():
            if node_id not in all_genomic_data:
                all_genomic_data[node_id] = regions
            else:
                all_genomic_data[node_id].extend(regions)
    return all_genomic_data

def is_overlap(region1, region2):
    """ Check if two regions overlap, considering chromosome IDs. """
    chrom1, start1, end1 = region1[0], region1[1], region1[2]
    chrom2, start2, end2 = region2[0], region2[1], region2[2]
    
    if chrom1 != chrom2:
        return False
    
    return max(start1, start2) < min(end1, end2)

def calculate_coverage_details(node_regions, node1, node2):
    """ Calculate detailed coverage information between two nodes. """
    total_regions_node1 = len(node_regions[node1])
    total_regions_node2 = len(node_regions[node2])
    overlapped_count = 0

    if total_regions_node1 == 0:
        return overlapped_count, total_regions_node1, total_regions_node2, 0

    for region1 in node_regions[node1]:
        for region2 in node_regions[node2]:
            if is_overlap(region1, region2):
                overlapped_count += 1
                break  # Stop checking once an overlap is found for this region

    coverage_percentage = (overlapped_count / total_regions_node1) * 100
    return overlapped_count, total_regions_node1, total_regions_node2, coverage_percentage

def calculate_coverage_and_output(genomic_data):
    """ Calculate and output detailed coverage for each node pair with size filtering. """
    nodes = sorted(genomic_data.keys(), key=lambda k: len(genomic_data[k]))
    
    for i, node1 in enumerate(nodes):
        for node2 in nodes[i+1:]:
            overlapped, total1, total2, coverage = calculate_coverage_details(genomic_data, node1, node2)
            print(f"{node1}\t{node2}\t{overlapped}\t{total1}\t{total2}\t{coverage:.2f}%")


file_paths = []
for filename in open("dum_ten_nodes.txt"):
    file_paths.append(filename.strip("\n"))
combined_genomic_data = merge_genomic_data(file_paths)
calculate_coverage_and_output(combined_genomic_data)

