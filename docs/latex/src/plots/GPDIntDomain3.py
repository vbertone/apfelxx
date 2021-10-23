import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

xi = 0.5

plt.axis("off")
plt.arrow(-1.2, -1, 2.3, 0, lw = 2, head_width = 0.05, head_length = 0.04, fc = "k")
plt.arrow(-1, -1.2, 0, 2.3, lw = 2, head_width = 0.03, head_length = 0.05, fc = "k")
plt.plot([0.7, 0.7], [-1.02, -0.98], lw = 1, c = "k")
plt.plot([-1.02, -0.98], [0.6, 0.6], lw = 1, c = "k")
plt.plot([-0.1, -0.1], [-1.02, -0.98], lw = 1, c = "k")

plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.text(-1.1, 1.1, r"$y$")
plt.text(1.1, -1.15, r"$x$")

plt.text(-1.1, 0.55, r"$1$")
plt.text(0.68, -1.15, r"$1$")
plt.text(-0.12, -1.15, r"$\xi$")

plt.text(0.1, 0.3, r"$\mathcal{P}_1$", color = "red", fontsize = 30)
plt.text(-0.83, 0, r"$\mathcal{P}_1$", color = "red", fontsize = 30)
plt.text(-0.65, 0, r"$+$", color = "black", fontsize = 30)
plt.text(-0.5, 0, r"$\mathcal{P}_2$", color = "blue", fontsize = 30)
plt.text(-0.65, 0.8, r"$\mathcal{P}_2$", color = "blue", fontsize = 30)
plt.text(0.1, 0.85, r"$\displaystyle y = \frac{1}{\xi}x$", fontsize = 18)

t  = np.linspace(0.9, 1.2, 100)
t2 = np.linspace(-0.15, 0.9, 100)
t3 = np.linspace(-1, -0.15, 100)
t4 = np.linspace(-1, 0.6, 100)

x,  y  = [-1, 0.7, -1, -1], [-1, 0.6, 0.6, -1]
x1,  y1  = [-1, -0.1, -0.1, -1], [-1, -0.15, 1.1, 1.1]

plt.plot([-0.1 for y in t], t, "k--", lw = 1.5)
plt.plot([-0.1 for y in t2], t2, "k-", lw = 1.5)

plt.plot([-0.1 for y in t3], t3, "k--", lw = 1)
plt.plot([0.7 for y in t4], t4, "k--", lw = 1)

ys = np.linspace(-1, 1.1, 100)
def xs(ys):
    return 0.9 * ( ys + 1 ) / 1.6 - 1

plt.plot(xs(ys), ys, "k-.", lw = 1)

plt.plot(x,  y,  c = "k", lw = 1.5)
plt.fill(x, y, alpha = 0.3, color = "red")
plt.fill(x1, y1, alpha = 0.3, color = "blue")

plt.savefig('GPDIntDomain3.pdf')
#plt.show()
