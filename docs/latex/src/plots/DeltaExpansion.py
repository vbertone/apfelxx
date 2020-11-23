"""
Plot demonstrating the integral as the area under a curve.

Although this is a simple example, it demonstrates some important tweaks:

    * A simple line plot with custom color and line width.
    * A shaded region created using a Polygon patch.
    * A text label with mathtext rendering.
    * figtext calls to label the x- and y-axes.
    * Use of axis spines to hide the top and right spines.
    * Custom tick placement and labels.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def func(x):
    return 5 / x

a, b = 10, 10  # integral limits
x = np.linspace(0, 10)
y = func(x)

fig, ax = plt.subplots()
plt.plot(x, y, 'r', linewidth=2)
plt.ylim(ymin=0)
plt.xlim(xmin=0)

# Make the shaded region
ix = np.linspace(a, b)
iy = func(ix)
verts = [(0, 0)] + [(a, 0)] + [(a, b)] + [(0, b)] 
poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
ax.add_patch(poly)

#plt.text(0.5 * (a + b), 30, r"$\int_a^b f(x)\mathrm{d}x$",
#         horizontalalignment='center', fontsize=20)

plt.figtext(0.9, 0.05, '$x$')
plt.figtext(0.1, 0.9, '$y$')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')

ax.set_xticks((a, b))
ax.set_xticklabels(('$a$', '$b$'))
ax.set_yticks([])

plt.savefig('DeltaExpansion.pdf')
plt.show()
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.patches import Polygon
import matplotlib.ticker as ticker

eps = 0.1
seps = 0.316227766016838

def func(x):
    return eps / x

def func1(x):
    return x

pl.style.use('seaborn-talk')
pl.rc('axes', linewidth=2)
ax = pl.gca()

#ax.get_yaxis().set_tick_params(which='both', direction='in', width=2, right='on')
#ax.get_xaxis().set_tick_params(which='both', direction='in', width=2, top='on')
#ax.get_yaxis().set_tick_params(which='major', size=10)
#ax.get_xaxis().set_tick_params(which='major', size=10)
#ax.get_yaxis().set_tick_params(which='minor', size=5)
#ax.get_xaxis().set_tick_params(which='minor', size=5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

pl.xlabel(r'$x$', fontsize=22, labelpad=-2, horizontalalignment='right', x=1.0)
pl.ylabel(r'$y$', fontsize=22, labelpad=-2, horizontalalignment='right', y=1.0)
pl.xlim(0, 1.2)
pl.ylim(0, 1.2)
pl.xticks(size=20)
pl.yticks(size=20)

ax.text(0.12, 1.1, r'$xy = \epsilon$', fontsize=22)
ax.text(0.95, 1.1,  r'$x = y$', fontsize=22)
ax.text(0.09, -0.08, r'$\epsilon$', fontsize=22)
ax.text(0.27, -0.08, r'$\sqrt{\epsilon}$', fontsize=22)

ax.text(-0.05, 0.09,  r'$\epsilon$', fontsize=22)
ax.text(-0.08, 0.305, r'$\sqrt{\epsilon}$', fontsize=22)

a, b = 1, 1  # integral limits

pl.axhline(y=eps, xmin=0, xmax=0.8333333333, linewidth=1, color = 'k', linestyle='--')
pl.axvline(x=eps, ymin=0, ymax=0.8333333333, linewidth=1, color = 'k', linestyle='--')

pl.axhline(y=seps, xmin=0, xmax=0.263523138347365, linewidth=1, color = 'k', linestyle='--')
pl.axvline(x=seps, ymin=0, ymax=0.263523138347365, linewidth=1, color = 'k', linestyle='--')

x = np.linspace(0, 1.2)
y = func(x)

pl.plot(x, y, 'r', linewidth=2)

x = np.linspace(0, 1.2)
y = func1(x)

pl.plot(x, y, 'b--', linewidth=2)

# Make the shaded region
ix = np.linspace(a, b)
iy = func(ix)
verts = [(0, 0)] + [(a, 0)] + [(a, b)] + [(0, b)] 
poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
ax.add_patch(poly)

#pl.plot(data[:,0], data[:,1], 'b-', label=r'Fixed order at $\mathcal{O}(\alpha_s)$')
#pl.plot(data[:,0], data[:,2], 'g-', label='Resummed at NLL')
#pl.plot(data[:,0], data[:,3], 'c-', label='Double counting')
#pl.plot(data[:,0], data[:,4], 'r-', label='Matched')

ax.get_xaxis().set_major_locator(ticker.MultipleLocator(1))
ax.get_yaxis().set_major_locator(ticker.MultipleLocator(1))

pl.legend(fontsize=22)

pl.savefig('DeltaExpansion.pdf')
pl.show()
