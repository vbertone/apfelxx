import ruamel.yaml as yaml
import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings
from scipy.interpolate import make_interp_spline, BSpline

# Loada data
data = np.loadtxt("FOvsAsy2.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [4, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set_title(r"\textbf{SIDIS at $\mathcal{O}(\alpha_s)$, $\sqrt{s}=10.5$ GeV}")
ax1.text(0.0002, 0.2, r"\textbf{$Q^2 = 2$ GeV$^2$}", fontsize = 16)
ax1.text(0.0002, 0.1, r"\textbf{$x = 0.1$}", fontsize = 16)
ax1.text(0.0002, 0.05, r"\textbf{$z = 0.2$}", fontsize = 16)
ax1.set(ylabel = r"$\displaystyle\left|\frac{d\sigma}{dy dz dQ dq_T}\right|$")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim([0.0001, 1])
ax1.set_ylim([0.0001, 10])
ax1.plot(data[:, 0], np.absolute(data[:, 1]), color = "red",  label = r"\textbf{Fixed order}")
ax1.plot(data[:, 0], np.absolute(data[:, 2]), color = "blue", label = r"\textbf{Asymptotic}")
ax1.plot(data[:, 0], np.absolute(data[:, 1] - data[:, 2]), color = "orange", label = r"\textbf{Difference}")
ax1.legend(fontsize = 20)

ax2.set_xlabel(r"\textbf{$q_T$ [GeV]}")
ax2.set_ylabel(r"\textbf{Ratio}", fontsize = 16)
ax2.set_ylim([0.55, 1.45])
ax2.plot(data[:, 0], np.absolute(data[:, 1] / data[:, 2]), color = "green")
ax2.plot(data[:, 0], np.absolute(data[:, 1] / data[:, 1]), color = "black", ls = "--", lw = 1.5)
ax2.set_xlim([0.0001, 1])

plt.savefig("FOvsAsy2.pdf")
plt.close()
