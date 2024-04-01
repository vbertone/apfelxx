import matplotlib.pyplot as plt
import MatplotlibSettings
import numpy as np
import scipy.special as sp
from scipy.integrate import quad

def psi(t):
    return t * np.tanh(np.pi * np.sinh(t) / 2)

def psip(t):
    return ( np.pi * t * np.cosh(t) + np.sinh(np.pi * np.sinh(t)) ) / ( 1 + np.cosh(np.pi * np.sinh(t)) )

# Parameters
qT = 1
h = min(0.005 * qT, 0.005)

def fb(b):
    return b * sp.j0(z) * np.exp(- b**2 / 2)
def fk(k):
    return np.exp(- k**2 / 2)

#a1 = 0.4
#a2 = 0.5
#def fb(b):
#    return b * sp.j0(b) * sp.j0(a1 * b) * sp.j0(a2 * b) * np.exp(- b**2 / 2)
#def fk(k):
#    return quad(fb, 0, 10, limit = 1000)[0]

#h = 0.01
#def fb(b):
#    return b / b * sp.j0(z)
#def fk(k):
#    return 1 / k

#h = 0.02
#def fb(b):
#    return b * sp.j0(z) / ( 1 + b**2 )
#def fk(k):
#    return sp.k0(k)

# Compute integral
nterm = 200
xi = sp.jn_zeros(0, nterm)
z  = np.pi * psi(h * xi / np.pi) / h
w  = ( np.pi / qT ) * ( sp.y0(xi) / sp.j1(xi) )
tm = w * psip(h * xi / np.pi) * fb(z / qT)
integ = np.array([sum(tm[:k]) for k in range(nterm)])

for k in range(nterm):
    print(k, integ[k], fk(qT), integ[k] / fk(qT))

plt.xlabel(r"$N$")
plt.ylabel(r"\textbf{Relative accuracy}")
plt.xlim(0, 30)
plt.ylim(0, 1.1)
plt.plot(range(nterm), integ / integ, c = "black", ls = "--", lw = 1.5)
plt.plot(range(nterm), integ / fk(qT), c = "blue", marker = "o", lw = 0)

plt.savefig("OgataTransformSingle.pdf")
plt.close()
