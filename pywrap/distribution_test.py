import apfelpy as ap
import numpy as np

# Class to define the analytical expression of LO splitting function
# P0qq
class p0qq(ap.Expression):
  def Regular(self, x):
    return - 2 * ap.constants.CF * ( 1 + x )
  def Singular(self, x):
    return 4 * ap.constants.CF / ( 1 - x )
  def Local(self, x):
    return 4 * ap.constants.CF * np.log( 1 - x ) + 3 * ap.constants.CF

# Class to define the analytical expression of the squared LO
# splitting function P0qq
class p0qq2(ap.Expression):
  def Regular(self, x):
    return 4 * ap.constants.CF * ap.constants.CF * ( - 4. * np.log(x) / ( 1. - x ) - 4. * ( 1. + x ) * np.log( 1. - x ) + 3. * ( 1. + x ) * np.log(x) - ( x + 5. ) )
  def Singular(self, x):
    return 4 * ap.constants.CF * ap.constants.CF * ( 8. * np.log( 1. - x ) + 6. ) / ( 1. - x )
  def Local(self, x):
    return 4 * ap.constants.CF * ap.constants.CF * ( 4. * np.log( 1. - x )**2 + 6 * np.log( 1. - x ) + 9. / 4. - 4. * np.pi**2 / 6. )

t = ap.Timer()

# Grid
g = ap.Grid([ap.SubGrid(80, 1e-5, 3), ap.SubGrid(50, 1e-1, 3), ap.SubGrid(40, 8e-1, 3)])

# Distribution
d = ap.Distribution(g, lambda x: 1 - x);

# Print distribution
d.Print()

# Expression
p = p0qq()

# Operator
print("\nInitialization ...")
t.start()
O = ap.Operator(g, p)
t.stop();

# Multiply operator by the distribution to create a new distribution
print("\nConvolution between operator and distribution (O * d) ...")
t.start();
Od = O * d;
t.stop();

# Multiply operator by itself to create a new operator
print("\nConvolution between two operators (O * O) ...")
t.start();
OO = O * O
t.stop();

# Check the numerical accuracy of "Od" by comparing with the analytical result
print("\nChecking the numerical accuracy of O * d ... ")
for ix in range(g.GetJointGrid().nx()):
      x = g.GetJointGrid().GetGrid()[ix]
      # Analytic result for x \int_x^1 dy Pqq(y) ( 1 - x / y )
      Ix = ap.constants.CF * ( - 2 * ( 3. / 2. - x - x**2 / 2. ) + 4 * ( 1 - x ) * np.log( 1 -  x ) + 3 * ( 1 - x ) + 2 * x * ( np.log(x) + 1 - x ) );
      print(format(x, ".6e"), "\t\t", format(Od.Evaluate(x), ".6e"), "\t\t", format(Ix, ".6e"), "\t\t", format(Od.Evaluate(x) / Ix, ".6e"))

# Check the numerical accuracy of "Od" by comparing with the
# analytical result Analytical expression of P0qq \otimes P0qq
p2 = p0qq2()
O2 = ap.Operator(g, p2)

# Now convolute both "OO" and "O2" with the test distribution "d" and
# check the result
OOd = OO * d
O2d = O2 * d

print("\nChecking the numerical accuracy of O * O ... ")
for ix in range(g.GetJointGrid().nx()):
    x = g.GetJointGrid().GetGrid()[ix]
    print(format(x, ".6e"), "\t\t", format(OOd.Evaluate(x), ".6e"), "\t\t", format(O2d.Evaluate(x), ".6e"), "\t\t", format(OOd.Evaluate(x) / O2d.Evaluate(x), ".6e"))
