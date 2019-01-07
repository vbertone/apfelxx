import numpy as np
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

data = np.loadtxt('FONLL.dat')

pl.title(r"SIDIS cross section", fontsize=22)

pl.xlabel(r'$q_T$ [GeV]', fontsize=22, labelpad=-2)
pl.ylabel(r'$\frac{d\sigma}{dx dy dz dq_T}$', fontsize=26, labelpad=-2.5)
pl.xlim(0.01, 100)
pl.ylim(0.01, 10)
pl.yscale('log')
pl.xscale('log')
pl.xticks(size=20)
pl.yticks(size=20)

ax.text(0.012, 0.2048, r'$\sqrt{S} = 300.9$ GeV', fontsize=22)
ax.text(0.012, 0.128,  r'$x = 10^{-3}$', fontsize=22)
ax.text(0.012, 0.08,   r'$ z = 10^{-1}$', fontsize=22)
ax.text(0.012, 0.05,  r'$y = 0.5$', fontsize=22)

ax.text(0.2, 0.02,  r'PDF set: CT14nlo', fontsize=22)
ax.text(0.2, 0.0125, r'FF set: DSS07 ($h+\overline{h}$)', fontsize=22)

pl.plot(data[:,0], data[:,1], 'b-', label=r'Fixed order at $\mathcal{O}(\alpha_s)$')
pl.plot(data[:,0], data[:,2], 'g-', label='Resummed at NLL')
pl.plot(data[:,0], data[:,3], 'c-', label='Double counting')
pl.plot(data[:,0], data[:,4], 'r-', label='Matched')

pl.legend(fontsize=22)

pl.savefig('FONLL.pdf')
pl.show()
