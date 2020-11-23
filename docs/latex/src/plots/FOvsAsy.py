import numpy as np
#import pylab as pl
import matplotlib.pyplot as pl

pl.style.use('seaborn-talk')
pl.rc('axes', linewidth=2)
ax = pl.gca()

ax.get_yaxis().set_tick_params(which='both', direction='in', width=2, right='on')
ax.get_xaxis().set_tick_params(which='both', direction='in', width=2, top='on')
ax.get_yaxis().set_tick_params(which='major', size=10)
ax.get_xaxis().set_tick_params(which='major', size=10)
ax.get_yaxis().set_tick_params(which='minor', size=5)
ax.get_xaxis().set_tick_params(which='minor', size=5)

data = np.loadtxt('FOvsAsy.dat')

pl.title(r"Fixed order vs. asymptotic limit at $\mathcal{O}(\alpha_s)$", fontsize=22)

pl.xlabel(r'$q_T$ [GeV]', fontsize=22, labelpad=-2)
pl.ylabel(r'$\frac{d\sigma}{dx dy dz dq_T}$', fontsize=26, labelpad=-2.5)
pl.xlim(0.01, 20)
pl.ylim(0.02, 10)
pl.yscale('log')
pl.xscale('log')
pl.xticks(size=20)
pl.yticks(size=20)

ax.text(0.012, 0.1024, r'$\sqrt{S} = 300.9$ GeV', fontsize=22)
ax.text(0.012, 0.064,  r'$x = 10^{-3}$', fontsize=22)
ax.text(0.012, 0.04,   r'$ z = 10^{-1}$', fontsize=22)
ax.text(0.012, 0.025,  r'$y = 0.5$', fontsize=22)

ax.text(0.7, 0.04,  r'PDF set: CT14nlo', fontsize=22)
ax.text(0.7, 0.025, r'FF set: DSS07 ($h+\overline{h}$)', fontsize=22)

pl.plot(data[:,0], data[:,1], 'r-', label='Fixed order')
pl.plot(data[:,0], data[:,2], 'b-', label='Asymptotic limit')

pl.legend(fontsize=22)

pl.savefig('FOvsAsy.pdf')
pl.show()
