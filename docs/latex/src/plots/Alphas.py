import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("Alphas.dat")

# Both
plt.xlabel(r"$\mu$ \textbf{[GeV]}")
plt.ylabel(r"$\alpha_s(\mu)$")
plt.xlim(0, 95)
plt.ylim(0.1, 0.4)
#plt.xscale("log")
plt.yscale("log")
plt.text(70, 0.199, r"\textbf{$\Downarrow$}")
plt.text(70, 0.251, r"\textbf{$\Downarrow$}")
plt.plot(data[:, 0], data[:, 1], c = "blue", linestyle = "--", label = r"\textbf{Analytic NLL solution}", lw = 1)
plt.plot(data[:, 0], data[:, 2], c = "blue", linestyle = "--", lw = 1)
plt.plot([91.2], [0.118],     "ro", markersize = 3, label = r"\textbf{1) $\alpha_s^{(n_f=5)}(M_Z)=0.118$}")
plt.plot([1], [0.3847609],    "yo", markersize = 3, label = r"\textbf{2) $\alpha_s^{(n_f=5)}(1\mbox{ GeV})=0.385$}")
plt.plot([91.2], [0.1164206], "go", markersize = 3, label = r"\textbf{3) $\alpha_s^{(n_f=5)}(M_Z)=0.116$}")

#plt.plot(data[:, 0], data[:, 3], c = "red", linestyle = "-", label = r"\textbf{Numeric}")
#plt.plot(data[:, 0], data[:, 4], c = "red", linestyle = "-")
lgnd = plt.legend(fontsize = 18, labelspacing = 2.5)

# Change the marker size manually for both lines
lgnd.legendHandles[1]._legmarker.set_markersize(10)
lgnd.legendHandles[2]._legmarker.set_markersize(10)
lgnd.legendHandles[3]._legmarker.set_markersize(10)

plt.savefig("Alphas.pdf")
plt.close()
