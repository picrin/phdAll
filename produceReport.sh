set -e
jupyter nbconvert masterBlood.ipynb --to latex
pdflatex masterBlood.tex
jupyter nbconvert normalise.ipynb --to latex
pdflatex normalise.tex
jupyter nbconvert analyseLimma.ipynb --to latex
pdflatex analyseLimma.tex
cd overleaf
xelatex main
biber main
xelatex main
mv main.pdf ../litReview.pdf
