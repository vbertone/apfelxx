#!/bin/bash

pandoc -s ../../README.md -o secs/README.rst
pandoc -s --toc ../latex/src/CCDIS.tex              -o secs/CCDIS.rst
pandoc -s --toc ../latex/src/CovarianceMatrix.tex   -o secs/CovarianceMatrix.rst
pandoc -s --toc ../latex/src/DIS.tex		    -o secs/DIS.rst
pandoc -s --toc ../latex/src/DrellYanTMD.tex	    -o secs/DrellYanTMD.rst
pandoc -s --toc ../latex/src/EvolDIS.tex	    -o secs/EvolDIS.rst
pandoc -s --toc ../latex/src/EvolutionCode.tex	    -o secs/EvolutionCode.rst
pandoc -s --toc ../latex/src/GPDs.tex		    -o secs/GPDs.rst
pandoc -s --toc ../latex/src/IntegralStucture.tex   -o secs/IntegralStucture.rst
pandoc -s --toc ../latex/src/MSbarMass.tex	    -o secs/MSbarMass.rst
pandoc -s --toc ../latex/src/MatchingConditions.tex -o secs/MatchingConditions.rst
pandoc -s --toc ../latex/src/SIDISTMD.tex	    -o secs/SIDISTMD.rst
pandoc -s --toc ../latex/src/TMDEvolution.tex	    -o secs/TMDEvolution.rst
pandoc -s --toc ../latex/src/Transversity.tex	    -o secs/Transversity.rst

cp ../latex/src/plots/*pdf secs/plots/.

make html
