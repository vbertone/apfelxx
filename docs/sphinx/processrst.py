import numpy as np
import glob

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

    # Close input file
    fin.close()

    # Write on the input file
    fout = open(s, "wt")
    fout.write(out)
    fout.close()
