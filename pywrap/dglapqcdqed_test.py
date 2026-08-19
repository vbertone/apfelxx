import apfelpy as ap
import numpy as np

# x-space grid
g = ap.Grid([ap.SubGrid(100, 1e-5, 3), ap.SubGrid(100, 1e-1, 3), ap.SubGrid(100, 6e-1, 3), ap.SubGrid(80, 8.5e-1, 5)])

# Initial scale
mu0 = np.sqrt(2)

# Final scale
mu = 100

# Vectors of thresholds
QuarkThresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]
LeptonThresholds = [0, 0, 1.777]

# Perturbative order
PerturbativeOrder = 2

# Running strong and electromagnetic couplings
a = ap.AlphaQCDQED(0.35, 7.496252e-3, np.sqrt(2), QuarkThresholds, LeptonThresholds, PerturbativeOrder)
Couplings = ap.TabulateObjectMatrix(a, 100, 0.9, 1001, 3)
Alphas  = lambda mu: Couplings.Evaluate(mu)(0, 0)
Alphaem = lambda mu: Couplings.Evaluate(mu)(1, 0)

# Initialize QCDxQED evolution objects
DglapObj   = ap.initializers.InitializeDglapObjectsQCDQED(g, QuarkThresholds, LeptonThresholds)
DglapObjOp = ap.initializers.InitializeDglapObjectsQCDQED(g, QuarkThresholds, LeptonThresholds, True)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFsQCDQED, mu0, PerturbativeOrder, Alphas, Alphaem)
EvolvedOps  = ap.builders.BuildDglap(DglapObjOp,                             mu0, PerturbativeOrder, Alphas, Alphaem)

# Tabulate PDFs and operators
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)
TabulatedOps  = ap.TabulateObjectSetO(EvolvedOps,  50, 1, 1000, 3)

# Evolve PDFs to the final scale
pdfs = ap.utilities.PlusMinusQCDQEDToPhys(EvolvedPDFs.Evaluate(mu).GetObjects())

# Interpolate the tabulated PDFs
tpdfs = ap.utilities.PlusMinusQCDQEDToPhys(TabulatedPDFs.Evaluate(mu).GetObjects())

# Set the appropriate convolution basis for the evolution operators and
# convolute them with the initial-scale distributions
tops = TabulatedOps.Evaluate(mu)
tops.SetMap(ap.EvolveDistributionsBasisQCDQED())
pdfs0 = ap.SetD(ap.EvolveDistributionsBasisQCDQED(), ap.utilities.DistributionMap(g, ap.utilities.LHToyPDFsQCDQED, mu0))
oppdfs = ap.utilities.PlusMinusQCDQEDToPhys((tops * pdfs0).GetObjects())

# Print results
xlha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1]

print("\nAlphaQCD(Q) = ", Alphas(mu))
print("AlphaQED(Q) = ", Alphaem(mu))

header = ("\n   x    ", "   u-ubar   ", "   d-dbar   ", " 2(ubr+dbr) ", "   c+cbar   ",
          "    gluon   ", "   photon   ", "   e^-+e^+  ", "  mu^-+mu^+ ", " tau^-+tau^+")


def line(d, x):
    return [(d[2] - d[-2]).Evaluate(x),
            (d[1] - d[-1]).Evaluate(x),
            2 * (d[-2] + d[-1]).Evaluate(x),
            (d[4] + d[-4]).Evaluate(x),
            d[0].Evaluate(x),
            d[22].Evaluate(x),
            2 * (d[11] + d[-11]).Evaluate(x),
            2 * (d[13] + d[-13]).Evaluate(x),
            2 * (d[15] + d[-15]).Evaluate(x)]



def scalar_line(d):
    return [d[2] - d[-2],
            d[1] - d[-1],
            2 * (d[-2] + d[-1]),
            d[4] + d[-4],
            d[0],
            d[22],
            2 * (d[11] + d[-11]),
            2 * (d[13] + d[-13]),
            2 * (d[15] + d[-15])]


def show(title, rows):
    print("\n" + title)
    print(*header)
    for x, vals in rows:
        print(format(x, ".1e"), " ", "  ".join(format(v, ".4e") for v in vals))


show("Direct Evolution:", [(x, line(pdfs, x)) for x in xlha])
show("Evolution through the evolution operator:", [(x, line(oppdfs, x)) for x in xlha])
show("Evolution through the interpolated evolution operator:",
     [(x, scalar_line(ap.utilities.PlusMinusQCDQEDToPhys((ap.SetD(ap.EvolveDistributionsBasisQCDQED(), tops.Evaluate(x).GetObjects()) * pdfs0).Squash()))) for x in xlha])
show("Interpolation on the PDF table (all x for each Q):", [(x, line(tpdfs, x)) for x in xlha])
show("Interpolation on the PDF table as a map (x and Q independently):",
     [(x, scalar_line(ap.utilities.PlusMinusQCDQEDToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu)))) for x in xlha])
