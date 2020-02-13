import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

x, y = [-1, 1, 0.6, 0.6, 1, -1, -0.6, -0.6, -1], [1, 1, 0.6, -0.6, -1, -1, -0.6, 0.6, 1]
x1, y1 = [-0.6, 0, -0.6], [0.6, 0, -0.6]
x2, y2 = [0.6, 0, 0.6], [0.6, 0, -0.6]

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

plt.text(0.64, 0.07, r"$\xi$")
plt.text(-0.76, 0.07, r"$-\xi$")

t = np.linspace(-1, 1, 100)
plt.plot(t, t, "k--", lw = 1.5)
plt.plot(t, -t, "k--", lw = 1.5)
plt.plot(x, y, c = "k", lw = 1.5)
plt.fill(x, y, alpha = 0.3)
plt.fill(x1, y1, alpha = 0.3, c = "r")
plt.fill(x2, y2, alpha = 0.3, c = "r")


plt.savefig('GPDIntDomain.pdf')
#plt.show()
