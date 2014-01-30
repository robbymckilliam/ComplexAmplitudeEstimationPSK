cd code/data/
mpost plotM2.mp
mpost plotM4.mp
mpost plotM8.mp
mpost plotncM2.mp
mpost plotncM4.mp
mpost plotncM8.mp
mpost plotturbo.mp
cd ../..

pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
