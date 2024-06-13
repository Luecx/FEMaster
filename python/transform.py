import argparse
from tqdm import tqdm

def process_input_deck(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    processed_lines = []
    temp_line = ""

    for line in tqdm(lines, desc="Processing lines"):
        entries = line.strip().split(',')
        if len(entries) >= 8:
            temp_line += ",".join(entries)
        else:
            processed_line = temp_line + ",".join(entries)
            processed_lines.append(processed_line.strip(",") + "\n")
            temp_line = ""

    if temp_line:
        processed_lines.append(temp_line.strip(","))

    with open(output_file, 'w') as file:
        file.writelines(processed_lines)
    print(f"Processing complete. Output written to {output_file}")

# Setting up argument parser
parser = argparse.ArgumentParser(description="Process an input deck file.")
parser.add_argument("input_file", type=str, help="Path to the input file")
parser.add_argument("output_file", type=str, help="Path to the output file")

# Parse arguments
args = parser.parse_args()

# Process the file
process_input_deck(args.input_file, args.output_file)

