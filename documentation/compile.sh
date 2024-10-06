#!/bin/bash

# Name of the LaTeX file (without extension)
file_name="document"

# Compile the LaTeX file using pdflatex (twice to resolve references)
pdflatex "$file_name.tex"
pdflatex "$file_name.tex"

# Remove auxiliary and temporary files generated during compilation
rm -f *.aux *.log *.out *.toc *.bbl *.blg
rm -f pages/*.aux

# Print success message
echo "Compilation complete. Temporary files deleted, and final PDF is ready: $file_name.pdf"
