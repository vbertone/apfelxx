import apfelpy as ap
import numpy as np

# x-space grid
g = ap.Grid([ap.SubGrid(100,1e-5,3), ap.SubGrid(60,1e-1,3), ap.SubGrid(50,6e-1,3), ap.SubGrid(50,8e-1,3)])

# Initial scale
mu0 = np.sqrt(2)

# Final scale
mu = 100

# Vectors of masses and thresholds
Thresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Perturbative order
PerturbativeOrder = 2

# Running coupling
a = ap.AlphaQCD(0.35, np.sqrt(2), Thresholds, PerturbativeOrder)
Alphas = ap.TabulateObject(a, 100, 0.9, 1001, 3)

# Initialize QCD evolution objects
DglapObj   = ap.initializers.InitializeDglapObjectsQCD(g, Thresholds)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)

# Compute results
tpdfs = ap.utilities.QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects())

# Print results
xlha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1]
print("\nAlphaQCD(Q) = " , Alphas.Evaluate(mu))
print("\n   x    "
          , "   u-ubar   "
          , "   d-dbar   "
          , " 2(ubr+dbr) "
          , "   c+cbar   "
          , "    gluon   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format((tpdfs[2] - tpdfs[-2]).Evaluate(x), ".4e"), " ",
              format((tpdfs[1] - tpdfs[-1]).Evaluate(x), ".4e"), " ",
              format(2 * (tpdfs[-2] + tpdfs[-1]).Evaluate(x), ".4e"), " ",
              format((tpdfs[4] + tpdfs[-4]).Evaluate(x), ".4e"), " ",
              format(tpdfs[0].Evaluate(x), ".4e"))
