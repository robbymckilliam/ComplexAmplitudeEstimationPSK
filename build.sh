cd code/data/
mpost plotM2.mp
mpost plotM4.mp
mpost plotM8.mp
mpost plotncM2.mp
mpost plotncM4.mp
mpost plotncM8.mp
cd ../..

pdflatex paper1.tex
pdflatex paper2.tex
bibtex paper1
bibtex paper2
pdflatex paper1.tex
pdflatex paper2.tex
pdflatex paper1.tex
pdflatex paper2.tex