@echo off
REM Name of the LaTeX file (without extension)
set file_name=document

REM Step 1: Compile the LaTeX file using pdflatex (first run to generate .aux files)
pdflatex %file_name%.tex

REM Step 2: Run biber to process the bibliography
biber %file_name%

REM Step 3: Compile the LaTeX file again to integrate bibliography
pdflatex %file_name%.tex

REM Step 4: Final pdflatex run to resolve cross-references
pdflatex %file_name%.tex

REM Remove auxiliary and temporary files generated during compilation
del /Q %file_name%.aux %file_name%.log %file_name%.out %file_name%.toc %file_name%.bbl %file_name%.blg %file_name%.bcf %file_name%.run.xml
del /Q pages\*.aux

REM Print success message
echo Compilation complete. Temporary files deleted, and final PDF is ready: %file_name%.pdf
pause
