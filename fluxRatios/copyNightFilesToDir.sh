#!/bin/bash

# Check if correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <csv_file> <output_directory>"
    exit 1
fi

# Input CSV file containing file paths
CSV_FILE=$1

# Output directory where files will be copied
OUTPUT_DIR=$2

# Check if the CSV file exists
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file does not exist!"
    exit 1
fi

# Check if the output directory exists, if not create it
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory does not exist. Creating it..."
    mkdir -p "$OUTPUT_DIR"
fi

# Read CSV file line by line and copy the files
while IFS=, read -r FILE_PATH; do
    # Strip any leading/trailing whitespaces (if necessary)
    FILE_PATH=$(echo "$FILE_PATH" | xargs)
    
    # Check if the file exists
    if [ -f "$FILE_PATH" ]; then
        echo "Copying $FILE_PATH to $OUTPUT_DIR"
        cp "$FILE_PATH" "$OUTPUT_DIR"
    else
        echo "File not found: $FILE_PATH"
    fi
done < "$CSV_FILE"

echo "File copying completed!"
