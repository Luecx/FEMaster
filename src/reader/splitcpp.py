import re
import os

# Input file path
input_file = "reader_process.cpp"

# Output directory for the split files
output_dir = "reader_parts"
os.makedirs(output_dir, exist_ok=True)

# Read the entire input file
with open(input_file, "r") as file:
    content = file.read()

# Regular expression to match each `Reader::process_` function
function_regex = re.compile(
    r"void\s+Reader::(process_\w+)\s*\((.*?)\)\s*\{", re.DOTALL
)

# Extract the namespace from the input file
namespace_regex = re.compile(r"namespace\s+fem::reader\s*\{", re.DOTALL)
namespace_match = namespace_regex.search(content)

if not namespace_match:
    print("Namespace 'fem::reader' not found.")
    exit(1)

namespace_declaration = "#include \"reader.h\"\n\nnamespace fem::reader {\n"
namespace_closing = "} // namespace fem::reader\n"

# Function to extract content for each function
matches = list(function_regex.finditer(content))

if not matches:
    print("No functions found.")
    exit(1)

# Iterate over matches to split the content
for i, match in enumerate(matches):
    function_name = match.group(1)
    start_index = match.start()
    end_index = matches[i + 1].start() if i + 1 < len(matches) else len(content)

    function_content = content[start_index:end_index].strip()

    # Create the new file name
    output_file = os.path.join(output_dir, f"reader_{function_name}.cpp")

    # Write the function content to the new file
    with open(output_file, "w") as out_file:
        out_file.write(namespace_declaration)
        out_file.write(function_content + "\n")
        out_file.write(namespace_closing)

print(f"Split {len(matches)} functions into separate files in '{output_dir}' directory.")