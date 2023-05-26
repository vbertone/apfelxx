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

# CKM matrix elements
def fCKM(Q):
    return ap.constants.CKM2

# Initialize QCD evolution objects
DglapObj = ap.initializers.InitializeDglapObjectsQCD(g, Thresholds)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)

# Initialize coefficient functions
F2PlusCCObj  = ap.initializers.InitializeF2CCPlusObjectsZM(g,  Thresholds)
F2MinusCCObj = ap.initializers.InitializeF2CCMinusObjectsZM(g, Thresholds)
FLPlusCCObj  = ap.initializers.InitializeFLCCPlusObjectsZM(g,  Thresholds)
FLMinusCCObj = ap.initializers.InitializeFLCCMinusObjectsZM(g, Thresholds)
F3PlusCCObj  = ap.initializers.InitializeF3CCPlusObjectsZM(g,  Thresholds)
F3MinusCCObj = ap.initializers.InitializeF3CCMinusObjectsZM(g, Thresholds)

# Initialize structure functions
F2p = ap.builders.BuildStructureFunctions(F2PlusCCObj,  TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)
F2m = ap.builders.BuildStructureFunctions(F2MinusCCObj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)
FLp = ap.builders.BuildStructureFunctions(FLPlusCCObj,  TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)
FLm = ap.builders.BuildStructureFunctions(FLMinusCCObj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)
F3p = ap.builders.BuildStructureFunctions(F3PlusCCObj,  TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)
F3m = ap.builders.BuildStructureFunctions(F3MinusCCObj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fCKM)

# Tabulate Structure functions
F2total  = ap.TabulateObjectD(lambda Q: F2p[0].Evaluate(Q) - F2m[0].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F2light  = ap.TabulateObjectD(lambda Q: F2p[1].Evaluate(Q) - F2m[1].Evaluate(Q) + F2p[2].Evaluate(Q) - F2m[2].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F2charm  = ap.TabulateObjectD(lambda Q: F2p[4].Evaluate(Q) - F2m[4].Evaluate(Q) + F2p[5].Evaluate(Q) - F2m[5].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F2bottom = ap.TabulateObjectD(lambda Q: F2p[3].Evaluate(Q) - F2m[3].Evaluate(Q) + F2p[6].Evaluate(Q) - F2m[6].Evaluate(Q), 50, 1, 200, 3, Thresholds)

FLtotal  = ap.TabulateObjectD(lambda Q: FLp[0].Evaluate(Q) - FLm[0].Evaluate(Q), 50, 1, 200, 3, Thresholds)
FLlight  = ap.TabulateObjectD(lambda Q: FLp[1].Evaluate(Q) - FLm[1].Evaluate(Q) + FLp[2].Evaluate(Q) - FLm[2].Evaluate(Q), 50, 1, 200, 3, Thresholds)
FLcharm  = ap.TabulateObjectD(lambda Q: FLp[4].Evaluate(Q) - FLm[4].Evaluate(Q) + FLp[5].Evaluate(Q) - FLm[5].Evaluate(Q), 50, 1, 200, 3, Thresholds)
FLbottom = ap.TabulateObjectD(lambda Q: FLp[3].Evaluate(Q) - FLm[3].Evaluate(Q) + FLp[6].Evaluate(Q) - FLm[6].Evaluate(Q), 50, 1, 200, 3, Thresholds)

F3total  = ap.TabulateObjectD(lambda Q: F3p[0].Evaluate(Q) - F3m[0].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F3light  = ap.TabulateObjectD(lambda Q: F3p[1].Evaluate(Q) - F3m[1].Evaluate(Q) + F3p[2].Evaluate(Q) - F3m[2].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F3charm  = ap.TabulateObjectD(lambda Q: F3p[4].Evaluate(Q) - F3m[4].Evaluate(Q) + F3p[5].Evaluate(Q) - F3m[5].Evaluate(Q), 50, 1, 200, 3, Thresholds)
F3bottom = ap.TabulateObjectD(lambda Q: F3p[3].Evaluate(Q) - F3m[3].Evaluate(Q) + F3p[6].Evaluate(Q) - F3m[6].Evaluate(Q), 50, 1, 200, 3, Thresholds)

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
              format(2 * F2light.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F2charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F2bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F2total.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  FLlight   "
          , "  FLcharm   "
          , "  FLbottom  "
          , "  FLtotal   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(2 * FLlight.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * FLcharm.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * FLbottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * FLtotal.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  F3light   "
          , "  F3charm   "
          , "  F3bottom  "
          , "  F3total   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(2 * F3light.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F3charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F3bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(2 * F3total.EvaluatexQ(x, Q), ".4e"))

t.start()

k = 1000000
print("\nInterpolating ", k, " times a single F2 on the grid... ")
t.start()
for i in range(k):
  F2total.EvaluatexQ(0.05, Q)
t.stop()
