import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

xi = 0.5

x, y = [-1, 1, xi, xi, 1, -1, -xi, -xi, -1], [1, 1, xi, -xi, -1, -1, -xi, xi, 1]
x1, y1 = [-xi, 0, -xi], [xi, 0, -xi]
x2, y2 = [xi, 0, xi], [xi, 0, -xi]
x3, y3 = [-xi, xi, xi, -xi], [1, 1, -1, -1]

plt.axis("off")
plt.arrow(-1.2, 0, 2.3, 0, lw = 2, head_width = 0.05, head_length = 0.04, fc = "k")
plt.arrow(0, -1.2, 0, 2.3, lw = 2, head_width = 0.03, head_length = 0.05, fc = "k")
plt.plot([1, 1], [-0.02, 0.02], lw = 1, c = "k")
plt.plot([-1, -1], [-0.02, 0.02], lw = 1, c = "k")

plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.text(0.05, 1.1, r"$x'$")
plt.text(1.1, -0.15, r"$x$")

plt.text(0.98, 0.07, r"$1$")
plt.text(-1.05, 0.07, r"$-1$")
plt.text(-0.06, 0.88, r"$1$")
plt.text(-0.15, -0.98, r"$-1$")

plt.text(xi + 0.02, 0.07, r"$|\xi|$")
plt.text(- xi - 0.22, 0.07, r"$-|\xi|$")

plt.text(xi+(1-xi-0.33)/2, xi+(1-xi)/2, r"$A$", color = "blue")
plt.text(-xi-(1-xi-0.13)/2, xi+(1-xi)/2, r"$\overline{A}$", color = "blue")

plt.text(xi+(1-xi-0.33)/2, -xi-(1-xi+0.2)/2, r"$\overline{A}$", color = "blue")
plt.text(-xi-(1-xi-0.13)/2, -xi-(1-xi+0.2)/2, r"$A$", color = "blue")

plt.text((xi-0.13)/2, xi+(1-xi-0.33)/2, r"$B$", color = "red")
plt.text(-(xi+0.07)/2, xi+(1-xi-0.33)/2, r"$\overline{B}$", color = "red")

plt.text((xi-0.13)/2, -xi-(1-xi-0.13)/2, r"$\overline{B}$", color = "red")
plt.text(-(xi+0.07)/2, -xi-(1-xi-0.13)/2, r"$B$", color = "red")

plt.text(xi/2, (xi - 0.33)/2, r"$C$", color = "red")
plt.text(-(xi+0.2)/2, (xi - 0.33)/2, r"$\overline{C}$", color = "red")

plt.text(-(xi+0.2)/2, -(xi - 0.13)/2, r"$C$", color = "red")
plt.text(xi/2, -(xi - 0.13)/2, r"$\overline{C}$", color = "red")

t = np.linspace(-1, 1, 100)
#plt.plot(t, t, "k--", lw = 1.5)
#plt.plot(t, -t, "k--", lw = 1.5)

plt.plot([xi for y in t], t, "k--", lw = 1.5)
plt.plot([-xi for y in t], t, "k--", lw = 1.5)

plt.plot(t, t, "k--", lw = 1.5)
plt.plot(t, -t, "k--", lw = 1.5)

plt.plot(x, y, c = "k", lw = 1.5)
plt.fill(x, y, alpha = 0.3)
#plt.fill(x1, y1, alpha = 0.3, c = "r")
#plt.fill(x2, y2, alpha = 0.3, c = "r")
plt.fill(x3, y3, alpha = 0.3, c = "r")

plt.savefig('GPDIntDomain.pdf')
#plt.show()
