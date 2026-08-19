import apfelpy as ap
import numpy as np

# x-space grid
g = ap.Grid([ap.SubGrid(100, 1e-5, 3), ap.SubGrid(100, 1e-1, 3), ap.SubGrid(100, 6e-1, 3), ap.SubGrid(80, 8.5e-1, 5)])

# Initial scale
mu0 = np.sqrt(2)

# Final scale
mu = 100

# Vector of thresholds
Thresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Perturbative order
PerturbativeOrder = 2

# Running strong coupling
a = ap.AlphaQCD(0.35, np.sqrt(2), Thresholds, PerturbativeOrder)
Alphas = ap.TabulateObject(a, 100, 0.9, 1001, 3)

# Initialize QCD evolution objects in the physical basis
DglapObj   = ap.initializers.InitializeDglapObjectsQCDPhys(g, Thresholds)
DglapObjOp = ap.initializers.InitializeDglapObjectsQCDPhys(g, Thresholds, True)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFsPlusMinus, mu0, PerturbativeOrder, Alphas.Evaluate)
EvolvedOps  = ap.builders.BuildDglap(DglapObjOp,                               mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs and operators
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)
TabulatedOps  = ap.TabulateObjectSetO(EvolvedOps,  50, 1, 1000, 3)

# Evolve PDFs to the final scale
pdfs = ap.utilities.PlusMinusToPhys(EvolvedPDFs.Evaluate(mu).GetObjects())

# Interpolate the tabulated PDFs
tpdfs = ap.utilities.PlusMinusToPhys(TabulatedPDFs.Evaluate(mu).GetObjects())

# Set the appropriate convolution basis for the evolution operators and
# convolute them with the initial-scale distributions
tops = TabulatedOps.Evaluate(mu)
tops.SetMap(ap.PhysicalEvolveDistributionsBasisQCD())
pdfs0 = ap.SetD(ap.PhysicalEvolveDistributionsBasisQCD(), ap.utilities.DistributionMap(g, ap.utilities.LHToyPDFsPlusMinus, mu0))
oppdfs = ap.utilities.PlusMinusToPhys((tops * pdfs0).GetObjects())

# Print results
xlha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1]

print("\nAlphaQCD(Q) = ", Alphas.Evaluate(mu))

header = ("\n   x    ", "   u-ubar   ", "   d-dbar   ", " 2(ubr+dbr) ", "   c+cbar   ", "    gluon   ")


def line(d, x):
    return [(d[2] - d[-2]).Evaluate(x),
            (d[1] - d[-1]).Evaluate(x),
            2 * (d[-2] + d[-1]).Evaluate(x),
            (d[4] + d[-4]).Evaluate(x),
            d[0].Evaluate(x)]


def scalar_line(d):
    return [d[2] - d[-2],
            d[1] - d[-1],
            2 * (d[-2] + d[-1]),
            d[4] + d[-4],
            d[0]]


def show(title, rows):
    print("\n" + title)
    print(*header)
    for x, vals in rows:
        print(format(x, ".1e"), " ", "  ".join(format(v, ".4e") for v in vals))


show("Direct Evolution:", [(x, line(pdfs, x)) for x in xlha])
show("Evolution through the evolution operator:", [(x, line(oppdfs, x)) for x in xlha])
show("Evolution through the interpolated evolution operator:",
     [(x, scalar_line(ap.utilities.PlusMinusToPhys((ap.SetD(ap.PhysicalEvolveDistributionsBasisQCD(), tops.Evaluate(x).GetObjects()) * pdfs0).Squash()))) for x in xlha])
show("Interpolation on the PDF table (all x for each Q):", [(x, line(tpdfs, x)) for x in xlha])
show("Interpolation on the PDF table as a map (x and Q independently):",
     [(x, scalar_line(ap.utilities.PlusMinusToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu)))) for x in xlha])
