cd code/data/
mpost plotncM2.mp
mpost plotncM4.mp
mpost plotncM8.mp
cd ../..

pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex