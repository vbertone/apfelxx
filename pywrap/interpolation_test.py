import apfelpy as ap

# Test distribution
def xg(x):
    return x * ( 1 - x )

# Derivative of the distribution
def dxg(x):
    return 1 - 2 * x

# Integral of the distribution
def ixg(x):
    return x * x * ( 1. / 2 - x / 3 )

# Grid
g = ap.Grid([ap.SubGrid(100, 9.9e-6, 3), ap.SubGrid(100, 1e-1, 3), ap.SubGrid(40, 8e-1, 5)])

# Interpolated distribution
xgluon = ap.Distribution(g, xg)

# Test values
x = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.51, 0.6, 0.7, 0.8, 0.9]

print("\n       x         "
          + "   an. func.         "
          + "    inter. func.      "
          + "        ratio           "
          + "      an. deriv.        "
          + "     inter. deriv.      "
          + "        ratio           "
          + "      an. integ.        "
          + "     inter. integ.      "
          + "        ratio           ")

for ix in x:
    original = xg(ix)
    derorig  = dxg(ix)
    intorig  = ixg(0.65) - ixg(ix)
    interpol = xgluon.Evaluate(ix)
    derive   = xgluon.Derive(ix)
    integr   = xgluon.Integrate(ix, 0.65)
    print(format(ix, ".8e"), "\t",
              format(original, ".8e"), "\t",
              format(interpol, ".8e"), "\t",
              format(interpol / original, ".8e"), "\t",
              format(derorig, ".8e"), "\t",
              format(derive, ".8e"), "\t",
              format(derive / derorig, ".8e"), "\t",
              format(intorig, ".8e"), "\t",
              format(integr, ".8e"), "\t",
              format(integr / intorig, ".8e"))
print("\n")

t = ap.Timer()
nint = 1000000;
print("Performance test (", nint, " interpolations)... ")
for i in range(nint):
    xgluon.Evaluate(0.1111)
t.stop()
t.start()
print("Performance test (", nint, " derivations)... ")
for i in range(nint):
    xgluon.Derive(0.1111)
t.stop()
t.start()
print("Performance test (", 1000, " integrations)... ")
for i in range(1000):
    xgluon.Integrate(0.1111, 0.55555);
t.stop()
