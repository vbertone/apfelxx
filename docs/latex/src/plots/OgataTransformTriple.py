import matplotlib.pyplot as plt
import MatplotlibSettings
import numpy as np
import scipy.special as sp
from scipy.integrate import quad

def Weights(nterm, a1, a2):
    # Short-hand notation for the Bessel functions
    def J0(x):
        return sp.j0(x)
    def J1(x):
        return sp.j1(x)
    def Y0(x):
        return sp.y0(np.absolute(x))
    def K0(x):
        return sp.k0(np.absolute(x))
    def I0(x):
        return sp.i0(x)

    # Get first nterm zeros
    xi  = sp.jn_zeros(0, nterm)

    # Vector of scaling factors
    a1eff = 1e-16 if a1 == 0 or a1 == 1 else a1
    a2eff = 1e-16 if a2 == 0 or a2 == 1 or a2 == a1 else a2
    av = np.array([1, a1eff, a2eff])

    # Compute sorted pairs of zero's and weighting factors
    xib = np.array([])
    an  = np.array([])
    kn  = np.array([])
    for i in range(len(av)):
        xib = np.append(xib, xi / av[i])
        an  = np.append(an, i * np.ones(nterm))
        kn  = np.append(kn, np.arange(nterm))

    # Sort vectors in ascending order of xib
    Xi  = np.array(sorted(list(zip(xib, an, kn))))
    xib = np.sort(xib)

    # Triple Bessel integral
    def L(k, c1, c2):
        if (c1 + c2 <= 1):
            return np.pi * Y0(xi[k]) * J0(c1 * xi[k]) * J0(c2 * xi[k])
        elif (np.absolute(c1 - c2) >= 1):
            return 0
        else:
            def wftTotrhs(z):
                return z * ( K0(c1 * z) * K0(c2 * z) * I0(z) + K0(z) * K0(c1 * z) * I0(c2 * z) + K0(c2 * z) * K0(z) * I0(c1 * z) ) / ( z**2 + xi[k]**2 )
            return np.pi / 2 * ( J0(c1 * xi[k]) * J0(c2 * xi[k]) + Y0(c1 * xi[k]) * Y0(c2 * xi[k]) ) * Y0(xi[k]) + 4 / np.pi**2 * quad(wftTotrhs, 0, 100, limit = 1000)[0]

    # Compute weights
    wi = np.array([])
    for i in range(len(Xi)):
        xibi = Xi[i][0]
        ni   = int(Xi[i][1])
        ki   = int(Xi[i][2])
        nip1 = (ni + 1) % 3
        nip2 = (ni + 2) % 3
        if av[ni] - np.absolute(av[nip1] - av[nip2]) <= 0:
            wi = np.append(wi, 0)
        else:
            wi = np.append(wi, L(ki, av[nip1] / av[ni], av[nip2] / av[ni]) / ( av[ni] * J1(av[ni] * xibi) * J0(av[nip1] * xibi) * J0(av[nip2] * xibi) ))

    # Return coordinates and weights
    return xib, wi

# Change of variable functions
def psi(tv):
    return np.array([t if t > 10 else t * np.tanh(np.pi * np.sinh(t) / 2) for t in tv])

def psip(tv):
    return np.array([1 if t > 10 else ( np.pi * t * np.cosh(t) + np.sinh(np.pi * np.sinh(t)) ) / ( 1 + np.cosh(np.pi * np.sinh(t)) ) for t in tv])

# Finally compute integral
hT0 = [1, 1,   1, 1,   1,   1,   1,  1,   1,   1]
hT1 = [0, 0.5, 1, 1.5, 0.3, 0.7, 1,  1.2, 0.3, 1.2]
hT2 = [0, 0,   0, 0,   0.3, 0.7, 1,  1.2, 0.8, 1.5]
clr = ["red", "blue", "orange", "green", "purple", "brown", "pink", "gray", "olive", "cyan"]

for i in range(len(hT0)):
    # Integrand
    def fb(b):
        return b * sp.j0(hT0[i] * b) * sp.j0(hT1[i] * b) * sp.j0(hT2[i] * b) * np.exp(- b**2 / 2)
    
    hTv = - np.sort(-np.array([hT0[i], hT1[i], hT2[i]]))
    h = 0.01 * hTv[0]
    nterm = 10

    xib, wi = Weights(nterm, hTv[1] / hTv[0], hTv[2] / hTv[0])
    z  = np.pi * psi(h * xib / np.pi) / h / hTv[0]
    tm = wi * psip(h * xib / np.pi) * fb(z) / hTv[0]
    integ = np.array([sum(tm[:k]) for k in range(len(tm))])
    bfint = quad(fb, 0, 10 * hTv[0], limit = 1000)[0]

    k = 0
    for numi in integ:
        print(k, numi, bfint, numi / bfint)
        k += 1

    plt.xlabel(r"$N$")
    plt.ylabel(r"\textbf{Relative accuracy}")
    plt.xlim(0, 30)
    plt.ylim(-0.1, 1.3)
    if i == 0:
        plt.plot(range(len(integ)), integ / integ, c = "black", ls = "--", lw = 1.5)
    plt.plot(range(len(integ)), integ / bfint, c = clr[i], marker = "o", lw = 0.5, ls = "-", label = r"\textbf{$h_{T0} = " + str(hT0[i]) + "$, $h_{T1} = " + str(hT1[i]) + r"$, $h_{T2} = " + str(hT2[i]) + r"$}")
    plt.legend(fontsize = 14)

plt.savefig("OgataTransformTriple.pdf")
plt.close()
