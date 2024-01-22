import sys

def parse_line(line):
    parts = line.split()
    node_id, chromosome_id, start, end, _ = parts
    start, end = int(start), int(end)
    # Ensure start is less than end
    if start > end:
        start, end = end, start
    return node_id, chromosome_id, start, end

def merge_regions(regions):
    regions.sort(key=lambda x: (x[0], x[1], x[2]))

    merged = []
    current = list(regions[0])

    for region in regions[1:]:
        if region[0] == current[0] and region[1] == current[1] and region[2] <= current[3]:
            current[3] = max(current[3], region[3])
        else:
            merged.append(tuple(current))
            current = list(region)

    merged.append(tuple(current))
    return merged

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.txt output_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    regions = []
    with open(input_file, 'r') as file:
        for line in file:
            regions.append(parse_line(line))

    merged_regions = merge_regions(regions)

    with open(output_file, 'w') as file:
        for region in merged_regions:
            node_id, chromosome_id, start, end = region
            length = end - start
            file.write(f"{node_id}\t{chromosome_id}\t{start}\t{end}\t{length}\n")

if __name__ == "__main__":
    main()

