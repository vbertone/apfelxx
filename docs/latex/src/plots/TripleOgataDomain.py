import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

xi = 0.5

plt.axis("off")
plt.arrow(-1.2, -1, 2.3, 0, lw = 2, head_width = 0.05, head_length = 0.04, fc = "k")
plt.arrow(-1, -1.2, 0, 2.3, lw = 2, head_width = 0.03, head_length = 0.05, fc = "k")

plt.plot([0, 0], [-1.02, -0.98], lw = 1, c = "k")
plt.plot([-1.02, -0.98], [0, 0], lw = 1, c = "k")

plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.text(-1.15, 1.1, r"$c_2$")
plt.text(1.1, -1.15, r"$c_1$")

plt.text(-1.1, -0.05, r"$1$")
plt.text(-0.02, -1.15, r"$1$")

#plt.text(0.1, 0.3, r"$\mathcal{P}_1$", color = "red", fontsize = 30)
#plt.text(-0.83, 0, r"$\mathcal{P}_1$", color = "red", fontsize = 30)
#plt.text(-0.65, 0, r"$+$", color = "black", fontsize = 30)
#plt.text(-0.5, 0, r"$\mathcal{P}_2$", color = "blue", fontsize = 30)
#plt.text(-0.65, 0.8, r"$\mathcal{P}_2$", color = "blue", fontsize = 30)

c1 = np.linspace(-1, 0, 100)
plt.plot(c1, - c1 - 1, "k-", lw = 1)

c1 = np.linspace(0, 0.8, 100)
plt.plot(c1, c1 - 1, "k-", lw = 1)

c1 = np.linspace(0, 1, 100)
plt.plot(c1, c1 - 1, "k--", lw = 1)

c1 = np.linspace(-1, -0.2, 100)
plt.plot(c1, c1 + 1, "k-", lw = 1)

c1 = np.linspace(-1, 0, 100)
plt.plot(c1, c1 + 1, "k--", lw = 1)

x,  y  = [-1, 0, -1], [0, 1, 1]
plt.fill(x, y, alpha = 0.3, color = "blue")

x,  y  = [0, 1, 1], [-1, -1, 0]
plt.fill(x, y, alpha = 0.3, color = "blue")

x,  y  = [-1, 0, -1], [-1, -1, 0]
plt.fill(x, y, alpha = 0.3, color = "red")

x,  y  = [-1, 0, 1, 1, 0], [0, -1, 0, 1, 1]
plt.fill(x, y, alpha = 0.3, color = "green")

plt.savefig("TripleOgataDomain.pdf")
#plt.show()
