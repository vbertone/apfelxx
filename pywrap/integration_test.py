import apfelpy as ap
import numpy as np

print("- 1D integration")

# Integrand function.
fGK = ap.Integrator(np.log, ap.IntegrationMethod.GAUSS_KRONROD)
fGL = ap.Integrator(np.log, ap.IntegrationMethod.GAUSS_LEGENDRE)

# Print true value.
print("True value: ", 2 * ( np.log(2) - 1 ))

# Integrate using Gauss-Legendre with a given accuracy.
res1 = fGL.integrate(0, 2, 1e-5)

# Integrate using Gauss-Legendre with a different accruracy.
res2 = fGL.integrate(0, 2, 1e-3)

# Print results.
print("Gauss-Legendre: ", format(res1, ".5f"), "  ", format(res2, ".5f"), "  ", format(res1 / res2, ".5f"))

# Integrate using Gauss-Kronrod with a given accuracy.
res6 = fGK.integrate(0, 2, 1e-5)

# Integrate using Gauss-Kronrod with a different accruracy.
res7 = fGK.integrate(0, 2, 1e-3)

# Print results.
print("Gauss-Kronrod:  ", format(res6, ".5f"), "  ", format(res7, ".5f"), "  ", format(res6 / res7, ".5f"))

# Now revert integration and try integrate with and without fixed
# points.
res4 = fGK.integrate(0, 2, 1e-5)
res5 = fGK.integrate(0, 2, [-1, 0.1, 0.5, 0.7, 1.1, 3.4], 1e-5)

# Print results.
print("W/o and w/ fixed points: ", format(res4, ".5f"), "  ", format(res5, ".5f"), "  ", format(res4 / res5, ".5f"))

# Performance test
k = 10000
print("Integrating ", k, " times with Gauss-Legendre... ")
t = ap.Timer()
for i in range(k):
    fGL.integrate(0, 2, 1e-12)
t.stop()

print("Integrating ", k, " times with Gauss-Kronrod... ")
t.start()
for i in range(k):
    fGK.integrate(0, 2, 1e-12)
t.stop()

# 2D integration
print("\n- 2D integration")

# Integrand function.
f2dGK = ap.Integrator2D(lambda x, y: np.log(x) * np.log(y), ap.IntegrationMethod.GAUSS_KRONROD)
f2dGL = ap.Integrator2D(lambda x, y: np.log(x) * np.log(y), ap.IntegrationMethod.GAUSS_LEGENDRE)

# Print true value.
print("True value: ", 4 * pow(np.log(2) - 1, 2))

# Integrate using Gauss-Legendre with a given accuracy.
res2d1 = f2dGL.integrate(0, 2, 0, 2, 1e-7)

# Integrate using Gauss-Legendre with a different accruracy.
res2d2 = f2dGL.integrate(0, 2, 0, 2, 1e-3)

# Print results.
print("Gauss-Legendre: ", format(res2d1, ".5f"), "  ", format(res2d2, ".5f"), "  ", format(res2d1 / res2d2, ".5f"))

# Integrate using Gauss-Kronrod with a given accuracy.
res2d3 = f2dGK.integrate(0, 2, 0, 2, 1e-7)

# Integrate using Gauss-Kronrod with a different accruracy.
res2d4 = f2dGK.integrate(0, 2, 0, 2, 1e-3)

# Print results.
print("Gauss-Kronrod:  ", format(res2d3, ".5f"), "  ", format(res2d4, ".5f"), "  ", format(res2d3 / res2d4, ".5f"))

# Integrate nested method
def In1(x):
  fin = ap.Integrator(lambda y: np.log(x) * np.log(y))
  return fin.integrate(0, 2, 1e-7)
fnest1 = ap.Integrator(In1)
res2d5 = fnest1.integrate(0, 2, 1e-7)

# Integrate nested method
def In2(x):
  fin = ap.Integrator(lambda y: np.log(x) * np.log(y))
  return fin.integrate(0, 2, 1e-3)
fnest2 = ap.Integrator(In2)
res2d6 = fnest1.integrate(0, 2, 1e-3)

# Print results.
print("Nested:         ", format(res2d5, ".5f"), "  ", format(res2d6, ".5f"), "  ", format(res2d5 / res2d6, ".5f"))

# Performance test
k = 10
print("Integrating ", k, " times with Gauss-Legendre... ")
t.start()
for i in range(k):
  f2dGL.integrate(0, 2, 0, 2, 1e-7)
t.stop()

print("Integrating ", k, " times with Gauss-Kronrod... ")
t.start()
for i in range(k):
  f2dGK.integrate(0, 2, 0, 2, 1e-7)
t.stop()

print("Integrating ", k, " times with nested function... ")
t.start()
for i in range(k):
  fnest1.integrate(0, 2, 1e-7)
t.stop()
