#!/bin/bash

# Introduction
pandoc -s ../../README.md -o secs/README.rst

# Declare an array of string with the document names
declare -a DocArray=("IntegralStucture" "Interpolation" "Statistics" "CCDIS" "Transversity" "TMDEvolution" "DrellYanTMD" "SIDISTMD" "DIS" "GPDs" "MatchingConditions" "MSbarMass" "EvolutionCode" "EvolDIS")

# Translate tex files into rst using pandoc
for doc in ${DocArray[@]}; do
    pandoc -s --toc --citeproc --bibliography=../latex/src/bibliography.bib ../latex/src/$doc.tex -o secs/$doc.rst
done

# Apply some corrections to the rst files
python3 processrst.py

# Produce the html page
make html
