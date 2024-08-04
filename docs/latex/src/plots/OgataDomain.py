import matplotlib.pyplot as plt
import MatplotlibSettings
import numpy as np

fig, ax = plt.subplots()
fig.patch.set_visible(False)
ax.axis("off")

plt.xlim(-10, 10)
plt.ylim(-2, 10)
plt.arrow(-9, 0, 18, 0, lw = 1.5, head_width = 0.4, head_length = 0.4, overhang = 0.4, color = "black")
plt.arrow(0, -5, 0, 14, lw = 1.5, head_width = 0.4, head_length = 0.5, overhang = 0.4, color = "black")
plt.text(9, 0.5, r"$\mbox{\textbf{Re}}(z)$")
plt.text(0.5, 9, r"$\mbox{\textbf{Im}}(z)$")
plt.text(9, 8.5, r"$\mathcal{C}$", bbox = dict(facecolor = "none", edgecolor = "black", pad = 10.0, linewidth = 1.5))

plt.scatter([5], [0], marker = "x", color = "black")
plt.text(4.5, -1, r"$\xi_k$")
plt.text(8, -0.8, r"$R$")
plt.text(-9, -0.8, r"$-R$")
plt.text(0.2, 6.5, r"$c$")
plt.text(8.8, 3.3, r"$C_1$")
plt.text(-1.3, 6.5, r"$C_2$")
plt.text(-9.7, 3.3, r"$C_3$")

plt.arrow(-8.5, 0.3, 8, 0, lw = 1, head_width = 0.2, head_length = 0.2, overhang = 0.4, color = "red")
plt.arrow(-0.5, 0.3, 4.7, 0, lw = 1, head_width = 0, head_length = 0, overhang = 0, color = "red")
plt.arrow(5.8, 0.3, 1.5, 0, lw = 1, head_width = 0.2, head_length = 0.2, overhang = 0.4, color = "red")
plt.arrow(7.3, 0.3, 1.2, 0, lw = 1, head_width = 0, head_length = 0, overhang = 0, color = "red")
plt.arrow(8.5, 0.3, 0, 3, lw = 1, head_width = 0.2, head_length = 0.2, overhang = 0.4, color = "red")
plt.arrow(8.5, 3.3, 0, 3, lw = 1, head_width = 0, head_length = 0, overhang = 0, color = "red")
plt.arrow(8.5, 6.3, -8, 0, lw = 1, head_width = 0.2, head_length = 0.2, overhang = 0.4, color = "red")
plt.arrow(0.5, 6.3, -9, 0, lw = 1, head_width = 0, head_length = 0, overhang = 0, color = "red")
plt.arrow(-8.5, 6.3, 0, -3, lw = 1, head_width = 0.2, head_length = 0.2, overhang = 0.4, color = "red")
plt.arrow(-8.5, 3.3, 0, -3, lw = 1, head_width = 0, head_length = 0, overhang = 0, color = "red")

# (x - x0)**2 + (y - y0)**2 == r**2
x0 = 5
y0 = 0.3
r = 0.8
x = np.linspace(x0 - r, x0 + r)
y = np.sqrt(r**2 - (x - x0)**2) + y0
ax.plot(x, y, ls = "-", lw = 1, color = "red")

plt.savefig("OgataDomain.pdf")
plt.close()
