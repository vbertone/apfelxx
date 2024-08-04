import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings
from scipy.special import j0

x = np.linspace(0, 50, 2000)

plt.xlabel(r"$x$")
plt.ylabel(r"$J_0(x)J_0(\sqrt{2}x)J_0(\sqrt{3}x)$")
plt.xlim(0, 50)
plt.ylim(-0.02, 0.02)
plt.plot(x, j0(x) * j0(np.sqrt(2) * x) * j0(np.sqrt(3) * x), c = "blue", linestyle = "-", lw = 2)

plt.plot(x, (2 / np.pi / x)**(3/2) / np.sqrt(np.sqrt(2) * np.sqrt(3)), c = "red", linestyle = "--", lw = 1.5)
plt.plot(x, - (2 / np.pi / x)**(3/2) / np.sqrt(np.sqrt(2) * np.sqrt(3)), c = "red", linestyle = "--", lw = 1.5)

plt.plot(x, x - x, c = "black", linestyle = "--", lw = 1.5)

plt.savefig("BesselProduct.pdf")
plt.close()
