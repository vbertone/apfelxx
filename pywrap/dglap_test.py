import apfelpy as ap
import numpy as np

# x-space grid
g = ap.Grid([ap.SubGrid(100,1e-5,3), ap.SubGrid(60,1e-1,3), ap.SubGrid(50,6e-1,3), ap.SubGrid(50,8e-1,3)])

# Initial scale
mu0 = np.sqrt(2)

# Vectors of masses and thresholds
Thresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Perturbative order
PerturbativeOrder = 2

# Running coupling
a = ap.AlphaQCD(0.35, np.sqrt(2), Thresholds, PerturbativeOrder)
Alphas = ap.TabulateObject(a, 100, 0.9, 1001, 3)

# Initialize QCD evolution objects
DglapObj   = ap.initializers.InitializeDglapObjectsQCD(g, Thresholds)
DglapObjOp = ap.initializers.InitializeDglapObjectsQCD(g, Thresholds, True)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj,   ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)
EvolvedOps  = ap.builders.BuildDglap(DglapObjOp,                         mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)

# Tabulate Operators
TabulatedOps = ap.TabulateObjectSetO(EvolvedOps, 50, 1, 1000, 3)

# Final scale
mu = 100

# Compute results
print("Direct evolution (4th order Runge-Kutta) from Q0 = " , mu0 , " GeV to Q = " , mu , " GeV... ")

# Evolve PDFs to the final Scale
t = ap.Timer()
pdfs = ap.utilities.QCDEvToPhys(EvolvedPDFs.Evaluate(mu).GetObjects())
t.stop()

print("Interpolation of the tabulated PDFs... ")
t.start()
tpdfs = ap.utilities.QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects())
t.stop()

print("Interpolation of the tabulated evolution operators... ")
t.start()
tops = TabulatedOps.Evaluate(mu)
t.stop()

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

print("Direct Evolution:")
for x in xlha:
    print(format(x, ".1e"), " ",
              format((pdfs[2] - pdfs[-2]).Evaluate(x), ".4e"), " ",
              format((pdfs[1] - pdfs[-1]).Evaluate(x), ".4e"), " ",
              format(2 * (pdfs[-2] + pdfs[-1]).Evaluate(x), ".4e"), " ",
              format((pdfs[4] + pdfs[-4]).Evaluate(x), ".4e"), " ",
              format(pdfs[0].Evaluate(x), ".4e"))

print("\nInterpolation on the PDF table (all x for each Q):")
for x in xlha:
    print(format(x, ".1e"), " ",
              format((tpdfs[2] - tpdfs[-2]).Evaluate(x), ".4e"), " ",
              format((tpdfs[1] - tpdfs[-1]).Evaluate(x), ".4e"), " ",
              format(2 * (tpdfs[-2] + tpdfs[-1]).Evaluate(x), ".4e"), " ",
              format((tpdfs[4] + tpdfs[-4]).Evaluate(x), ".4e"), " ",
              format(tpdfs[0].Evaluate(x), ".4e"))

print("\nInterpolation on the PDF table as a map (x and Q independently):")
for x in xlha:
    DistMap = ap.utilities.QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x,mu))
    print(format(x, ".1e"), " ",
              format(DistMap[2] - DistMap[-2], ".4e"), " ",
              format(DistMap[1] - DistMap[-1], ".4e"), " ",
              format(2 * (DistMap[-2] + DistMap[-1]), ".4e"), " ",
              format(DistMap[4] + DistMap[-4], ".4e"), " ",
              format(DistMap[0], ".4e"))

k = 1000000
print("\nInterpolating ", k, " times a single PDF on the (x,Q) grid... ")
t.start()
for i in range(k):
  TabulatedPDFs.EvaluatexQ(0, 0.05, mu)
t.stop()

k = 100000
print("Interpolating ", k, " times a map of PDFs on the (x,Q) grid... ")
t.start()
for i in range(k):
  TabulatedPDFs.EvaluateMapxQ(0.05, mu)
t.stop()
