import apfelpy as ap
import numpy as np

# Constructor of QGrid with type double
qg = ap.QGrid(50, 1, 1000, 3, [0, 0, 0, np.sqrt(2), 4.5, 175.])
qg.Print()

# AlphaQCD running coupling
als = ap.AlphaQCD(0.35, np.sqrt(2), [0, 0, 0, np.sqrt(2), 4.5, 175], 2)

# Tabulate AlphaQCD on a QGrid
gas = ap.TabulateObject(als, 50, 1, 1000, 3)

print("\nAccuracy test:")
nQ   = 20
Qmin = 1.1
Qmax = 999.
Step = np.exp( np.log( Qmax / Qmin ) / ( nQ - 1 ) )
Q = Qmin
print("Q         \t\tDirect            \t\tInterpolated          \t\tRatio")
for iQ in range(nQ):
    print(format(Q, ".8e"), "\t\t", format(als.Evaluate(Q), ".8e"), "\t\t", format(gas.Evaluate(Q), ".8e"), "\t\t", format(als.Evaluate(Q) / gas.Evaluate(Q), ".8e"))
    Q *=Step

print("\nPerformance test:")
t = ap.Timer()
nQ   = 100000
Step = np.exp( np.log( Qmax / Qmin ) / ( nQ - 1 ) )

t.start()
print("Direct calculation of ", nQ, " points... ")
Q = Qmin
for iQ in range(nQ):
    als.Evaluate(Q)
    Q *= Step
t.stop()

t.start()
print("Interpolated calculation of ", nQ, " points... ")
Q = Qmin
for iQ in range(nQ):
    gas.Evaluate(Q)
    Q *= Step
t.stop()
