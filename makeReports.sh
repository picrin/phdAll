set -xe
cd report
pdflatex boxplots
biber boxplots
pdflatex boxplots
pdflatex boxplots
