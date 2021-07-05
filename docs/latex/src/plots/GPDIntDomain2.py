import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

xi = 0.5

plt.axis("off")
plt.arrow(-1.2, 0, 2.3, 0, lw = 2, head_width = 0.05, head_length = 0.04, fc = "k")
plt.arrow(0, -1.2, 0, 2.3, lw = 2, head_width = 0.03, head_length = 0.05, fc = "k")
plt.plot([1, 1], [-0.02, 0.02], lw = 1, c = "k")
plt.plot([-1, -1], [-0.02, 0.02], lw = 1, c = "k")

plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.text(0.05, 1.1, r"$y$")
plt.text(1.1, -0.15, r"$x$")

plt.text(0.98, 0.07, r"$1$")
plt.text(-1.05, 0.07, r"$-1$")
plt.text(-0.06, 0.88, r"$1$")
plt.text(-0.15, -0.98, r"$-1$")

plt.text(xi + 0.02, 0.07, r"$|\xi|$")
plt.text(- xi - 0.22, 0.07, r"$-|\xi|$")

plt.text(xi+(1-xi-0.4)/2, 1.12, r"$P_{\rm NS}\left(\frac{x}{y},\frac{2\xi}{y}\right)$", color = "blue")
plt.text(-1, -1.2, r"$-P_{\rm NS}\left(\frac{x}{y},\frac{2\xi}{y}\right)$", color = "blue")
plt.text(-1, xi - 0.2, r"$P_{\rm NS}'\left(\frac{x}{y},\frac{2\xi}{y}\right)$", color = "red")
plt.text(xi + 0.01, -xi + 0.1, r"$-P_{\rm NS}'\left(\frac{x}{y},-\frac{2\xi}{y}\right)$", color = "orange")

t = np.linspace(-1, 1, 100)

x,  y  = [-1, 1, xi, xi, 1, -1, -xi, -xi, -1], [1, 1, xi, -xi, -1, -1, -xi, xi, 1]
x1, y1 = [xi, xi, 1], [xi, 1, 1]
x2, y2 = [-xi, -xi, -1], [-xi, -1, -1]
x3, y3 = [-xi, -xi, xi, xi], [-xi, 1, 1, xi]
x4, y4 = [-xi, -xi, xi, xi], [-xi, -1, -1, xi]

plt.plot([xi for y in t], t, "k--", lw = 1.5)
plt.plot([-xi for y in t], t, "k--", lw = 1.5)

plt.plot(t, t, "k--", lw = 1.5)
plt.plot(t, -t, "k--", lw = 1.5)

plt.plot(x,  y,  c = "k", lw = 1.5)
plt.fill(x1, y1, alpha = 0.3, color = "blue")
plt.fill(x2, y2, alpha = 0.3, color = "blue")
plt.fill(x3, y3, alpha = 0.3, color = "red")
plt.fill(x4, y4, alpha = 0.3, color = "orange")

plt.savefig('GPDIntDomain2.pdf')
#plt.show()
