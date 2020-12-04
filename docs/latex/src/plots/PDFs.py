import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

pdfs = np.loadtxt("PDFs.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.text(0.01, 180, r"\textbf{NLL evolution}", fontsize = 16)
ax1.text(0.01, 60, r"\textbf{MMHT2014nlo68cl}", fontsize = 16)
ax1.text(0.01, 20,  r"\textbf{$\mu_0=1$ GeV, $\mu=100$ GeV}", fontsize = 16)
#ax1.text(3, 0.8, r"\textbf{$\sqrt{s}=13$ TeV}")
#ax1.text(3, 0.4, r"\textbf{$Q = M_Z$, $y=0$}")
ax1.set(ylabel = r"$xg(x,\mu)$")
ax1.set_xlim([0.00001, 1])
ax1.set_ylim([0.000005, 1000])
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.plot(pdfs[:,0], pdfs[:,1], "r-", label = r"\textbf{DGLAP + $\alpha_s$ numerical}")
ax1.plot(pdfs[:,0], pdfs[:,2], "b-", label = r"\textbf{DGLAP + $\alpha_s$ analytical}")
ax1.plot(pdfs[:,0], pdfs[:,3], "k--", lw = 1, label = r"\textbf{LHAPDF}")
ax1.legend(fontsize = 18, loc = "lower left")

ax2.set(xlabel = r"\textbf{$x$}")
ax2.set_ylabel(r"\textbf{Ratio to num.}", fontsize = 16)
ax1.set_xlim([0.00001, 1])
ax2.set_ylim([0.8, 1.2])
ax2.set_xscale("log")
ax2.plot(pdfs[:,0], pdfs[:,1]/pdfs[:,1], "r--")
ax2.plot(pdfs[:,0], pdfs[:,2]/pdfs[:,1], "b--")
ax2.plot(pdfs[:,0], pdfs[:,3]/pdfs[:,1], "k--", lw = 1)

plt.savefig("PDFs.pdf")
plt.close()
