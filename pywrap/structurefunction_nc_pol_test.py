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

# Initialize QCD polarised evolution objects
DglapObj = ap.initializers.InitializeDglapObjectsQCDpol(g, Thresholds)

# Construct the DGLAP objects
EvolvedPDFs = ap.builders.BuildDglap(DglapObj, ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate PDFs
TabulatedPDFs = ap.TabulateObjectSetD(EvolvedPDFs, 50, 1, 1000, 3)

# Initialize coefficient functions
g4Obj = ap.initializers.Initializeg4NCObjectsZM(g, Thresholds)
gLObj = ap.initializers.InitializegLNCObjectsZM(g, Thresholds)
g1Obj = ap.initializers.Initializeg1NCObjectsZM(g, Thresholds)

# Initialize structure functions
g4 = ap.builders.BuildStructureFunctions(g4Obj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fDq)
gL = ap.builders.BuildStructureFunctions(gLObj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fDq)
g1 = ap.builders.BuildStructureFunctions(g1Obj, TabulatedPDFs.EvaluateMapxQ, PerturbativeOrder, Alphas.Evaluate, fBq)

# Tabulate Structure functions
g4total  = ap.TabulateObjectD(g4[0].Evaluate, 50, 1, 200, 3, Thresholds)
g4light  = ap.TabulateObjectD(lambda Q: g4[1].Evaluate(Q) + g4[2].Evaluate(Q) + g4[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
g4charm  = ap.TabulateObjectD(g4[4].Evaluate, 50, 1, 200, 3, Thresholds)
g4bottom = ap.TabulateObjectD(g4[5].Evaluate, 50, 1, 200, 3, Thresholds)

gLtotal  = ap.TabulateObjectD(gL[0].Evaluate, 50, 1, 200, 3, Thresholds)
gLlight  = ap.TabulateObjectD(lambda Q: gL[1].Evaluate(Q) + gL[2].Evaluate(Q) + gL[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
gLcharm  = ap.TabulateObjectD(gL[4].Evaluate, 50, 1, 200, 3, Thresholds)
gLbottom = ap.TabulateObjectD(gL[5].Evaluate, 50, 1, 200, 3, Thresholds)

g1total  = ap.TabulateObjectD(g1[0].Evaluate, 50, 1, 200, 3, Thresholds)
g1light  = ap.TabulateObjectD(lambda Q: g1[1].Evaluate(Q) + g1[2].Evaluate(Q) + g1[3].Evaluate(Q), 50, 1, 200, 3, Thresholds)
g1charm  = ap.TabulateObjectD(g1[4].Evaluate, 50, 1, 200, 3, Thresholds)
g1bottom = ap.TabulateObjectD(g1[5].Evaluate, 50, 1, 200, 3, Thresholds)

# Timer
t = ap.Timer()

# Final scale
Q = 100

# Print results
xlha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1]
print("\nAlphaQCD(Q) = " , Alphas.Evaluate(Q))

print("\n   x    "
          , "  g4light   "
          , "  g4charm   "
          , "  g4bottom  "
          , "  g4total   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(g4light.EvaluatexQ(x, Q), ".4e"), " ",
              format(g4charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(g4bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(g4total.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  gLlight   "
          , "  gLcharm   "
          , "  gLbottom  "
          , "  gLtotal   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(gLlight.EvaluatexQ(x, Q), ".4e"), " ",
              format(gLcharm.EvaluatexQ(x, Q), ".4e"), " ",
              format(gLbottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(gLtotal.EvaluatexQ(x, Q), ".4e"))

print("\n   x    "
          , "  g1light   "
          , "  g1charm   "
          , "  g1bottom  "
          , "  g1total   "
         )
for x in xlha:
    print(format(x, ".1e"), " ",
              format(g1light.EvaluatexQ(x, Q), ".4e"), " ",
              format(g1charm.EvaluatexQ(x, Q), ".4e"), " ",
              format(g1bottom.EvaluatexQ(x, Q), ".4e"), " ",
              format(g1total.EvaluatexQ(x, Q), ".4e"))

t.start()

k = 1000000
print("\nInterpolating ", k, " times a single g1 on the grid... ")
t.start()
for i in range(k):
  g1total.EvaluatexQ(0.05, Q)
t.stop()
