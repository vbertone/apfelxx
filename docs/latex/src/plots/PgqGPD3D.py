import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Skewness
xi = 0.5

# Define splitting function
@np.vectorize
def Pgq(y, x):
    k = xi / x
    CF = 4. / 3.
    k2 = k * k
    y2 = y * y
    P1 = 2 * CF * ( 1 + ( 1 - y )**2 - k2 * y2 ) / y / ( 1 - k2 * y2 )
    P2 = - 2 * CF * ( 1 - k )**2 / k / ( 1 - k2 * y2 )
    return ( ( P1 if y < 1 else 0 ) + ( P2 if k > 1 else 0 ) if y > x else 0 )

# Make data.
#X = np.arange(0.001, 1, 0.01)
#Y = np.arange(0.001, 2, 0.01)
#X, Y = np.meshgrid(X, Y)
#Z = Pgq(Y, X)

x = np.logspace(-4, 0, 100)
print(x)

plt.yscale("log")
plt.xscale("log")
plt.xlim(0.0001, 1)
plt.plot(x, Pgq(0.001, x),  c = "k", lw = 1)
plt.plot(x, Pgq(0.01, x),  c = "k", lw = 1)
plt.plot(x, Pgq(0.1, x),  c = "k", lw = 1)
plt.plot(x, Pgq(0.2, x),  c = "k", lw = 1)
plt.plot(x, Pgq(0.5, x),  c = "k", lw = 1)
plt.plot(x, Pgq(0.8, x),  c = "k", lw = 1)

plt.savefig('PgqGPD3D.pdf')
#plt.show()
