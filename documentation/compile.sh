#!/bin/bash

# Name of the LaTeX file (without extension)
file_name="document"

# Step 1: Compile the LaTeX file using pdflatex (first run to generate .aux files)
pdflatex "$file_name.tex"

# Step 2: Run biber to process the bibliography
biber "$file_name"

# Step 3: Compile the LaTeX file again to integrate bibliography
pdflatex "$file_name.tex"

# Step 4: Final pdflatex run to resolve cross-references
pdflatex "$file_name.tex"

# Remove auxiliary and temporary files generated during compilation
rm -f *.aux *.log *.out *.toc *.bbl *.blg *.bcf *.run.xml
rm -f pages/*.aux

# Print success message
echo "Compilation complete. Temporary files deleted, and final PDF is ready: $file_name.pdf"
