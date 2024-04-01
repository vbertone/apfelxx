import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import MatplotlibSettings
import numpy as np
import scipy.special as sp
from scipy.integrate import quad

# Bessel functions
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

# Parameters
nterm = 200
h = 0.001
xi = sp.jn_zeros(0, nterm)

# Test interpolation
print("Interpolation test")
x = 10

# Single Bessel
w = np.array([- J0(x) / J1(xik) * ( 2 * xik / ( x**2 - xik**2 ) ) for xik in xi])
fk = np.exp(- xi**2 / 100)
print("Single Bessel: true value = ", np.exp(- x**2 / 100), ", interpolated value = ", np.dot(w, fk))

# Triple Bessel
a1 = 20
a2 = 0.1
w0  = np.array([- J0(x) * J0(a1 * x) * J0(a2 * x) / J1(xik) / J0(a1 * xik) / J0(a2 * xik) * ( 2 * xik / ( x**2 - xik**2 ) ) for xik in xi])
fk0 = np.exp(- xi**2 / 100)
w1  = np.array([- J0(x) * J0(a1 * x) * J0(a2 * x) / a1 / J0(xik / a1) / J1(xik) / J0(a2 * xik / a1) * ( 2 * xik / a1 / ( x**2 - ( xik / a1 )**2 ) ) for xik in xi])
fk1 = np.exp(- (xi / a1)**2 / 100)
w2  = np.array([- J0(x) * J0(a1 * x) * J0(a2 * x) / a2 / J0(xik / a2) / J0(a1 * xik / a2) / J1(xik) * ( 2 * xik / a2 / ( x**2 - ( xik / a2 )**2 ) ) for xik in xi])
fk2 = np.exp(- (xi / a2)**2 / 100)
print("Triple Bessel: true value = ", np.exp(- x**2 / 100), ", interpolated value = ", np.dot(w0, fk0) + np.dot(w1, fk1) + np.dot(w2, fk2))

# Test integrals of the interpolating functions
print("\nWeight integration test")
k = 0

# Single Bessel
def wfs(x):
    return np.absolute(x) * J0(x) / ( x - xi[k] )
print("Single Bessel: numerical value = ", quad(wfs, -10000, 10000, limit = 10000)[0], ", analytic value = ", - np.pi * np.absolute(xi[k]) * Y0(xi[k]))

# Test values of c1 and c2
c1 = 1.2
c2 = 2

def wft1(x):
    return np.absolute(x) * ( J0(c1 * x) * J0(c2 * x) - Y0(c1 * x) * Y0(c2 * x) ) * J0(x) / ( x - xi[k] )
def wft1rhs(z):
    return z * K0(c1 * z) * K0(c2 * z) * I0(z) / ( z**2 + xi[k]**2 )
print("Triple Bessel (comb. 1): numerical value = ", "{:e}".format(quad(wft1, -1000, -0.00001, limit = 100000)[0] + quad(wft1, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wft1, xi[k] + 0.00001, 1000, limit = 100000)[0]), ", analytic value = ", "{:e}".format(- 8 * np.absolute(xi[k]) / np.pi**2 * quad(wft1rhs, 0, 100, limit = 1000)[0]))

def wft2(x):
    return np.absolute(x) * ( J0(x) * J0(c1 * x) - Y0(x) * Y0(c1 * x) ) * J0(c2 * x) / ( x - xi[k] )
def wft2rhs(z):
    return z * K0(z) * K0(c1 * z) * I0(c2 * z) / ( z**2 + xi[k]**2 )
print("Triple Bessel (comb. 2): numerical value = ", "{:e}".format(quad(wft2, -1000, -0.00001, limit = 100000)[0] + quad(wft2, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wft2, xi[k] + 0.00001, 1000, limit = 100000)[0]), ", analytic value = ", "{:e}".format(- np.pi * np.absolute(xi[k]) * Y0(xi[k]) * J0(c1 * xi[k]) * J0(c2 * xi[k]) - 8 * np.absolute(xi[k]) / np.pi**2 * quad(wft2rhs, 0, 100, limit = 1000)[0]))

def wft3(x):
    return np.absolute(x) * ( J0(c2 * x) * J0(x) - Y0(c2 * x) * Y0(x) ) * J0(c1 * x) / ( x - xi[k] )
def wft3rhs(z):
    return z * K0(c2 * z) * K0(z) * I0(c1 * z) / ( z**2 + xi[k]**2 )
print("Triple Bessel (comb. 3): numerical value = ", "{:e}".format(quad(wft3, -1000, -0.00001, limit = 100000)[0] + quad(wft3, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wft3, xi[k] + 0.00001, 1000, limit = 100000)[0]), ", analytic value = ", "{:e}".format(- np.pi * np.absolute(xi[k]) * J0(c1 * xi[k]) * Y0(xi[k]) * J0(c2 * xi[k])- 8 * np.absolute(xi[k]) / np.pi**2 * quad(wft3rhs, 0, 100, limit = 1000)[0]))

def wft4(x):
    return np.absolute(x) * ( J0(c1 * x) * J0(c2 * x) * J0(x) - Y0(c1 * x) * Y0(c2 * x) * J0(x) - Y0(x) * Y0(c1 * x) * J0(c2 * x) - Y0(c2 * x) * Y0(x) * J0(c1 * x) ) / ( x - xi[k] )
print("Triple Bessel (comb. 4): numerical value = ", "{:e}".format(quad(wft4, -1000, -0.00001, limit = 100000)[0] + quad(wft4, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wft4, xi[k] + 0.00001, 1000, limit = 100000)[0]), ", analytic value = ", "{:e}".format(- np.pi * np.absolute(xi[k]) * ( J0(c1 * xi[k]) * J0(c2 * xi[k]) - Y0(c1 * xi[k]) * Y0(c2 * xi[k]) ) * Y0(xi[k])))

def wftTot(x):
    return np.absolute(x) * J0(c1 * x) * J0(c2 * x) * J0(x) / ( x - xi[k] )
def wftTotrhs(z):
    return z * ( K0(c1 * z) * K0(c2 * z) * I0(z) + K0(z) * K0(c1 * z) * I0(c2 * z) + K0(c2 * z) * K0(z) * I0(c1 * z) ) / ( z**2 + xi[k]**2 )
print("Triple Bessel (total): numerical value = ", "{:e}".format(quad(wftTot, -1000, -0.00001, limit = 100000)[0] + quad(wftTot, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wftTot, xi[k] + 0.00001, 1000, limit = 100000)[0]), ", analytic value = ", "{:e}".format(- np.pi * np.absolute(xi[k]) / 2 * ( J0(c1 * xi[k]) * J0(c2 * xi[k]) + Y0(c1 * xi[k]) * Y0(c2 * xi[k]) ) * Y0(xi[k]) - 4 * np.absolute(xi[k]) / np.pi**2 * quad(wftTotrhs, 0, 100, limit = 1000)[0]))

# Total triple Bessel weight
def WTriple(c1, c2):
    if (c1 + c2 < 1):
        return np.pi * Y0(xi[k]) * J0(c1 * xi[k]) * J0(c2 * xi[k])
    elif (np.absolute(c1 - c2) > 1):
        return 0
    else:
        def wftTotrhs(z):
            return z * ( K0(c1 * z) * K0(c2 * z) * I0(z) + K0(z) * K0(c1 * z) * I0(c2 * z) + K0(c2 * z) * K0(z) * I0(c1 * z) ) / ( z**2 + xi[k]**2 )
        return np.pi / 2 * ( J0(c1 * xi[k]) * J0(c2 * xi[k]) + Y0(c1 * xi[k]) * Y0(c2 * xi[k]) ) * Y0(xi[k]) + 4 / np.pi**2 * quad(wftTotrhs, 0, 100, limit = 1000)[0]

print("\nTesting the full weights")
c1 = 0.5
c2 = 0
for i in range(20):
    def wft(x):
        return - np.absolute(x) * J0(c1 * x) * J0(c2 * x) * J0(x) / ( x - xi[k] ) / xi[k]
    nu = quad(wft, -1000, -0.00001, limit = 100000)[0] + quad(wft, 0.00001, xi[k] - 0.00001, limit = 100000)[0] + quad(wft, xi[k] + 0.00001, 1000, limit = 100000)[0]
    print("{:e}".format(c1), "{:e}".format(c2), "{:e}".format(nu), "{:e}".format(WTriple(c1, c2)))
    #c1 += 0.1
    c2 += 0.1

# Make 3D plot
c1v = np.arange(0.0001, 3, 0.1)
c2v = np.arange(0.0001, 3, 0.1)
w = np.array([])
for c2 in c2v:
    for c1 in c1v:
        w = np.append(w, WTriple(c1, c2))
w = np.reshape(w, (-1, len(c2v)))
c1v, c2v = np.meshgrid(c1v, c2v)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.view_init(20, 25)

# Plot the surface.
surf = ax.plot_surface(c1v, c2v, w, cmap = cm.coolwarm, linewidth = 0, antialiased = False)

# Customize the z axis.
ax.set_zlim(-0.1, 1.8)
ax.set_xlim(0, 3)
ax.set_ylim(0, 3)
ax.set_zlabel(r"\textbf{$\mathcal{L}_1(c_1,c_2)$}", labelpad = 15)
ax.set_xlabel(r"\textbf{$c_1$}", labelpad = 15)
ax.set_ylabel(r"\textbf{$c_2$}", labelpad = 15)
ax.set_box_aspect(aspect = None, zoom = 0.8)

#plt.show()
plt.savefig("Weight3D.pdf", bbox_inches = "tight", pad_inches = 0)
