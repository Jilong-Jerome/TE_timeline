#!/bin/bash

# Get file path from command line argument
file=$1

# Check if file name was provided
if [ -z "$file" ]; then
    echo "Usage: $0 [file_path]"
    exit 1
fi

# Check if file exists
if [ ! -f "$file" ]; then
    echo "File not found: $file"
    exit 1
fi

# Initialize a line number counter
line_number=0

# Process each line
while IFS=$'\t' read -r node sequence; do
    # Increment line number
    ((line_number++))

    # Skip the first line (header)
    if [ $line_number -eq 1 ]; then
        continue
    fi

    # Skip empty lines
    if [ -z "$node" ] || [ -z "$sequence" ]; then
        continue
    fi

    # Extract species from the sequence name (first field of the '|' delimited string)
    species=$(echo "$sequence" | cut -d '|' -f1)

    # Create species directory if it doesn't exist
    species_dir="./${species}"
    if [ ! -d "$species_dir" ]; then
        mkdir -p "$species_dir"
    fi

    # Append sequence to the node-specific file
    echo "$sequence" >> "${species_dir}/${node}_seq.txt"
done < "$file"

