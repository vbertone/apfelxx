import ruamel.yaml as yaml
import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("NewPrescription.dat")

# Both
plt.xlabel(r"$q_T$ \textbf{[GeV]}")
plt.ylabel(r"$\displaystyle \frac{d\sigma}{dQdydq_T}$ \textbf{[pb/GeV$^2$]}")
plt.text(2.5, 1.5, r"\textbf{N$^3$LL}")
plt.text(2.5, 1.0, r"\textbf{NNPDF3.1 NNLO luxQED}")
plt.text(2.5, 0.5, r"$Q=M_Z\,,\;y=0$")
plt.xlim(0, 20)
plt.ylim(0, 4)
plt.plot(data[0:50, 2], data[0:50, 3], c = "blue", linestyle = "--", dashes=(5, 1), label = r"\textbf{CSS}")
plt.plot(data[0:50, 2], data[0:50, 6], c = "red", linestyle = "--", dashes=(5, 1), label = r"\textbf{Alternative}")
plt.legend()
plt.savefig("NewPrescription_qT.pdf")
plt.close()

# CSS
plt.xlabel(r"$q_T$ \textbf{[GeV]}")
plt.ylabel(r"$\frac{d\sigma}{dQdydq_T}$")
plt.xlim(0, 20)
plt.ylim(0, 4)
plt.plot(data[0:50, 2], data[0:50, 3], c = "blue", linestyle = "--", dashes=(5, 1), label = r"$b_0/b_{\rm max} = 1$")
plt.plot(data[0:50, 2], data[0:50, 4], c = "blue", linestyle = "--", dashes=(4, 2), label = r"$b_0/b_{\rm max} = 1.5$")
plt.plot(data[0:50, 2], data[0:50, 5], c = "blue", linestyle = "--", dashes=(3, 3), label = r"$b_0/b_{\rm max} = 0.67$")
plt.legend()
plt.savefig("NewPrescription_qT_CSS.pdf")
plt.close()

# Alternative
plt.xlabel(r"$q_T$ \textbf{[GeV]}")
plt.ylabel(r"$\frac{d\sigma}{dQdydq_T}$")
plt.xlim(0, 20)
plt.ylim(0, 4)
plt.plot(data[0:50, 2], data[0:50, 6], c = "red", linestyle = "--", dashes=(5, 1), label = r"$b_0/b_{\rm max} = 1$")
plt.plot(data[0:50, 2], data[0:50, 7], c = "red", linestyle = "--", dashes=(4, 2), label = r"$b_0/b_{\rm max} = 1.5$")
plt.plot(data[0:50, 2], data[0:50, 8], c = "red", linestyle = "--", dashes=(3, 3), label = r"$b_0/b_{\rm max} = 0.67$")
plt.legend()
plt.savefig("NewPrescription_qT_noCSS.pdf")
plt.close()

# CSS in bT
plt.xlabel(r"$b_T$ \textbf{[GeV$^{-1}$]}")
plt.ylabel(r"$b_T\widetilde{\sigma}(b_T)$")
plt.xlim(0, 1.5)
plt.ylim(0, 4)
plt.plot(data[50:150, 2], data[50:150, 3], c = "blue", linestyle = "--", dashes=(5, 1), label = r"$b_0/b_{\rm max} = 1$")
plt.plot(data[50:150, 2], data[50:150, 4], c = "blue", linestyle = "--", dashes=(4, 2), label = r"$b_0/b_{\rm max} = 1.5$")
plt.plot(data[50:150, 2], data[50:150, 5], c = "blue", linestyle = "--", dashes=(3, 3), label = r"$b_0/b_{\rm max} = 0.67$")
plt.legend()
plt.savefig("NewPrescription_bT_CSS.pdf")
plt.close()

# Alternative in bT
plt.xlabel(r"$b_T$ \textbf{[GeV$^{-1}$]}")
plt.ylabel(r"$b_T\widetilde{\sigma}(b_T)$")
plt.xlim(0, 1.5)
plt.ylim(0, 4)
plt.plot(data[50:150, 2], data[50:150, 6], c = "red", linestyle = "--", dashes=(5, 1), label = r"$b_0/b_{\rm max} = 1$")
plt.plot(data[50:150, 2], data[50:150, 7], c = "red", linestyle = "--", dashes=(4, 2), label = r"$b_0/b_{\rm max} = 1.5$")
plt.plot(data[50:150, 2], data[50:150, 8], c = "red", linestyle = "--", dashes=(3, 3), label = r"$b_0/b_{\rm max} = 0.67$")
plt.legend()
plt.savefig("NewPrescription_bT_noCSS.pdf")
plt.close()


