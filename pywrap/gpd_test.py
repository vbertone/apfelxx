import apfelpy as ap
import numpy as np

# x-space grid
g = ap.Grid([ap.SubGrid(100,1e-5,3), ap.SubGrid(200,1e-1,3), ap.SubGrid(100,9e-1,3)])
#g = ap.Grid([ap.SubGrid(30,1e-5,3), ap.SubGrid(20,1e-1,3), ap.SubGrid(10,9e-1,3)])

# Initial scale
mu0 = np.sqrt(2)

# Vectors of masses and thresholds
Thresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Perturbative order
PerturbativeOrder = 0

# Skewness
xi = 0.3

# Running coupling
a = ap.AlphaQCD(0.35, np.sqrt(2), Thresholds, PerturbativeOrder)
Alphas = ap.TabulateObject(a, 100, 0.9, 1001, 3)

# Initialize QCD evolution objects
GpdObj   = ap.initializers.InitializeGpdObjects(g, Thresholds, xi)
#GpdObjOp = ap.initializers.InitializeGpdObjects(g, Thresholds, xi, True)

# Construct the DGLAP objects
EvolvedGPDs = ap.builders.BuildDglap(GpdObj,   ap.utilities.LHToyPDFs, mu0, PerturbativeOrder, Alphas.Evaluate)
#EvolvedOps  = ap.builders.BuildDglap(GpdObjOp,                         mu0, PerturbativeOrder, Alphas.Evaluate)

# Tabulate GPDs
TabulatedGPDs = ap.TabulateObjectSetD(EvolvedGPDs, 30, 1, 100, 3)

# Tabulate Operators
#TabulatedOps = ap.TabulateObjectSetO(EvolvedOps, 30, 1, 100, 3)

# Final scale
mu = 10

# Compute results
print("Direct evolution (4th order Runge-Kutta) from Q0 = " , mu0 , " GeV to Q = " , mu , " GeV... ")

# Evolve GPDs to the final Scale
t = ap.Timer()
gpds = ap.utilities.QCDEvToPhys(EvolvedGPDs.Evaluate(mu).GetObjects())
t.stop()

print("Interpolation of the tabulated GPDs... ")
t.start()
tgpds = ap.utilities.QCDEvToPhys(TabulatedGPDs.Evaluate(mu).GetObjects())
t.stop()

#print("Interpolation of the tabulated evolution operators... ")
#t.start()
#tops = TabulatedOps.Evaluate(mu)
#t.stop()

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

#print("Direct Evolution:")
#for x in xlha:
#    print(format(x, ".1e"), " ",
#              format((gpds[2] - gpds[-2]).Evaluate(x), ".4e"), " ",
#              format((gpds[1] - gpds[-1]).Evaluate(x), ".4e"), " ",
#              format(2 * (gpds[-2] + gpds[-1]).Evaluate(x), ".4e"), " ",
#              format((gpds[4] + gpds[-4]).Evaluate(x), ".4e"), " ",
#              format(gpds[0].Evaluate(x), ".4e"))

print("\nInterpolation on the GPD table (all x for each Q):")
for x in xlha:
    print(format(x, ".1e"), " ",
              format((tgpds[2] - tgpds[-2]).Evaluate(x), ".4e"), " ",
              format((tgpds[1] - tgpds[-1]).Evaluate(x), ".4e"), " ",
              format(2 * (tgpds[-2] + tgpds[-1]).Evaluate(x), ".4e"), " ",
              format((tgpds[4] + tgpds[-4]).Evaluate(x), ".4e"), " ",
              format(tgpds[0].Evaluate(x), ".4e"))

print("\nInterpolation on the GPD table as a map (x and Q independently):")
for x in xlha:
    DistMap = ap.utilities.QCDEvToPhys(TabulatedGPDs.EvaluateMapxQ(x,mu))
    print(format(x, ".1e"), " ",
              format(DistMap[2] - DistMap[-2], ".4e"), " ",
              format(DistMap[1] - DistMap[-1], ".4e"), " ",
              format(2 * (DistMap[-2] + DistMap[-1]), ".4e"), " ",
              format(DistMap[4] + DistMap[-4], ".4e"), " ",
              format(DistMap[0], ".4e"))

t.stop()
