import numpy as np
import glob
import re

# Get section files
secs = glob.glob("./secs/*rst")

# Run over section files
for s in secs:
    fin = open(s, "rt")

    #for each line in the input file
    out = ""
    labels = []
    for line in fin:
        if line.find("\label{") >= 0:
            out += line.replace("\label{", ":label: ").replace("}", "")
            # Collect labels
            labels.append(line[line.find("{") + 1: line.find("}")])
        else:
            out += line

    # Replace lables
    for l in labels:
        out = out.replace("(`[" + l + "] <#" + l + ">`__)", ":eq:`" + l + "`")

    # Remove blanck line after .. math::
    out = out.replace(".. math::\n\n", ".. math::\n")

    # Replace "\bm" with "\vec" because it is not rendered
    out = out.replace("\\bm", "\\vec")

    # Replace the folder "plots" to that pointing to the original files
    out = out.replace("plots/", "../../latex/src/plots/")

    # Place the word "References" where the list of references starts
    out = out.replace(".. container:: references csl-bib-body hanging-indent", "**References**\n\n.. container:: references csl-bib-body hanging-indent")

    # Comment out :alt: to make figured render
    out = out.replace("   :alt:", ".. :alt:")

    # Close input file
    fin.close()

    # Write on the input file
    fout = open(s, "wt")
    fout.write(out)
    fout.close()
