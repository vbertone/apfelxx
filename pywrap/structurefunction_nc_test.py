import apfelpy as ap
import numpy as np

# Show banner
ap.Banner()

# x-space grid
g = ap.Grid([ap.SubGrid(100, 1e-5, 3), ap.SubGrid(60, 1e-1, 3), ap.SubGrid(50, 6e-1, 3), ap.SubGrid(50, 8e-1, 3)])

# Initial scale
mu0 = np.sqrt(2)

# Vectors of masses and thresholds
Thresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Perturbative order
PerturbativeOrder = 2

# Running coupling
a = ap.AlphaQCD(0.35, np.sqrt(2), Thresholds, PerturbativeOrder)
Alphas = ap.TabulateObject(a, 100, 0.9, 1001, 3)

# Effective charges
def fBq(Q):
    return ap.utilities.ElectroWeakCharges(Q, False)

def fDq(Q):
    return ap.utilities.ParityViolatingElectroWeakCharges(Q, False)

# Initialize QCD evolution objects
DglapObj = ap.initializers.InitializeDglapObjectsQCD(g, Thresholds)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)

# Initialize coefficient functions
F2Obj = ap.initializers.InitializeF2NCObjectsZM(g, Thresholds)
FLObj = ap.initializers.InitializeFLNCObjectsZM(g, Thresholds)
F3Obj = ap.initializers.InitializeF3NCObjectsZM(g, Thresholds)

# Initialize structure functions
F2 = ap.builders.BuildStructureFunctions(F2Obj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fBq)
FL = ap.builders.BuildStructureFunctions(FLObj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fBq)
F3 = ap.builders.BuildStructureFunctions(F3Obj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fDq)

# Tabulate Structure functions
F2total  = ap.TabulateObjectD(F2[0].Evaluate, 50, 1, 200, 3, Thresholds)
F2light  = ap.TabulateObjectD(lambda Q: F2[1].Evaluate(Q) + F2[2].Evaluate(Q) + F2[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F2charm  = ap.TabulateObjectD(F2[4].Evaluate, 50, 1, 200, 3, Thresholds)
F2bottom = ap.TabulateObjectD(F2[5].Evaluate, 50, 1, 200, 3, Thresholds)

FLtotal  = ap.TabulateObjectD(FL[0].Evaluate, 50, 1, 200, 3, Thresholds)
FLlight  = ap.TabulateObjectD(lambda Q: FL[1].Evaluate(Q) + FL[2].Evaluate(Q) + FL[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
FLcharm  = ap.TabulateObjectD(FL[4].Evaluate, 50, 1, 200, 3, Thresholds)
FLbottom = ap.TabulateObjectD(FL[5].Evaluate, 50, 1, 200, 3, Thresholds)

F3total  = ap.TabulateObjectD(F3[0].Evaluate, 50, 1, 200, 3, Thresholds)
F3light  = ap.TabulateObjectD(lambda Q: F3[1].Evaluate(Q) + F3[2].Evaluate(Q) + F3[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F3charm  = ap.TabulateObjectD(F3[4].Evaluate, 50, 1, 200, 3, Thresholds)
F3bottom = ap.TabulateObjectD(F3[5].Evaluate, 50, 1, 200, 3, Thresholds)

# Timer
t = ap.Timer()

# Final scale
Q = 100

# Print results
xlha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1]
print("\nAlphaQCD(Q) = " , Alphas.Evaluate(Q))

print("\n   x    "
          , "  F2light   "
          , "  F2charm   "
          , "  F2bottom  "
          , "  F2total   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(F2light.EvaluatexQ(x, Q), ".4e"), " ",
              format(F2charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(F2bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(F2total.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  FLlight   "
          , "  FLcharm   "
          , "  FLbottom  "
          , "  FLtotal   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(FLlight.EvaluatexQ(x, Q), ".4e"), " ",
              format(FLcharm.EvaluatexQ(x, Q), ".4e"), " ",
              format(FLbottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(FLtotal.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  F3light   "
          , "  F3charm   "
          , "  F3bottom  "
          , "  F3total   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(F3light.EvaluatexQ(x, Q), ".4e"), " ",
              format(F3charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(F3bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(F3total.EvaluatexQ(x, Q), ".4e"))

t.start()

k = 1000000
print("\nInterpolating ", k, " times a single F2 on the grid... ")
t.start()
for i in range(k):
  F2total.EvaluatexQ(0.05, Q)
t.stop()
