#!/bin/bash

# Declare an array of string with the document names
declare -a DocArray=("IntegralStucture" "Interpolation" "Statistics" "StructureFunctions" "CCDIS" "Transversity" "TMDEvolution" "DrellYanTMD" "SIDISTMD"  "GPDs" "MatchingConditions" "MSbarMass" "EvolutionCode" "EvolDIS" "SIDISCollinear" "QCD_QED_common_basis")

# Produce pdfs
for doc in ${DocArray[@]}; do
    pdflatex $doc
    bibtex   $doc
    pdflatex $doc
    pdflatex $doc
done

# Move pdfs and clean up
mv ./*.pdf ../pdf/.
rm ./*.aux ./*.log ./*.out ./*.toc ./*.idx ./*.blg ./*.bbl
