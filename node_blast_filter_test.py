import sys

def merge_overlapping_hits(hits):
    sorted_hits = sorted(hits, key=lambda x: x[1])
    merged_hits = []

    current_hit = sorted_hits[0]
    for hit in sorted_hits[1:]:
        if hit[1] <= current_hit[2]:
            # Merge hits
            current_hit = (current_hit[0], min(current_hit[1], hit[1]), max(current_hit[2], hit[2]), current_hit[3])
        else:
            # Add the completed merged hit to the list
            merged_hits.append(current_hit)
            current_hit = hit
    merged_hits.append(current_hit)  # Add the last hit

    return merged_hits

def filter_and_merge_blast_hits(input_file_path, output_file_path):
    hits_by_subject = {}
    
    with open(input_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.split()
            subject_id = fields[7]
            subject_start = min(int(fields[3]),int(fields[4]))
            subject_end = max(int(fields[3]),int(fields[4]))
            alignment_length = int(fields[6])
            query_id = fields[0]
            query_length = abs(int(query_id.split('|')[-1]))
            node = query_id.split('|')[0]  # Extract the node information

            if alignment_length >= 0.8 * query_length:
                if subject_id not in hits_by_subject:
                    hits_by_subject[subject_id] = []
                hits_by_subject[subject_id].append((subject_id, subject_start, subject_end, node))

    merged_hits = []
    for subject_id, hits in hits_by_subject.items():
        merged_hits.extend(merge_overlapping_hits(hits))

    with open(output_file_path, 'w') as output_file:
        for hit in merged_hits:
            merged_region_length = hit[2] - hit[1] + 1  # Calculate merged region length
            reformatted_hit = f"{hit[3]}\t{hit[0]}\t{hit[1]}\t{hit[2]}\t{merged_region_length}\n"
            output_file.write(reformatted_hit)

def main():
    if len(sys.argv) != 3:
        print("Usage: python blast_filter.py <input_blast_results.txt> <output_filtered_results.txt>")
        print("This script filters, merges overlapping hits, and reformats BLAST results, including node information and merged region length.")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_and_merge_blast_hits(input_file, output_file)

if __name__ == "__main__":
    main()

