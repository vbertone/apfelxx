#!/bin/bash

for ((index=0; index<2; index++))
do
    pdflatex manual.tex
    pdflatex IntegralStucture.tex
done

mv *.pdf ../pdf/.
rm *.aux *.log *.out *.toc *.idx
