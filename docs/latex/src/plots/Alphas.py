import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("Alphas.dat")

# Both
plt.xlabel(r"$\mu$ \textbf{[GeV]}")
plt.ylabel(r"$\alpha_s(\mu)$")
#plt.text(2.5, 1.5, r"\textbf{NLL}")
plt.xlim(0.8, 100)
plt.ylim(0.1, 0.4)
plt.xscale("log")
plt.yscale("log")
plt.plot(data[:, 0], data[:, 1], c = "blue", linestyle = "-", label = r"\textbf{Analytic NLL solution}", lw = 1.5)
plt.plot(data[:, 0], data[:, 2], c = "blue", linestyle = "-", lw = 1.5)
#plt.plot(data[:, 0], data[:, 3], c = "red", linestyle = "-", label = r"\textbf{Numeric}")
#plt.plot(data[:, 0], data[:, 4], c = "red", linestyle = "-")
plt.legend()
plt.savefig("Alphas.pdf")
plt.close()
