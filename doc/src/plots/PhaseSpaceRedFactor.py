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

data = np.loadtxt('PhaseSpaceRedFactor.dat')

pl.title(r"$Q = M_Z$, lepton cuts $p_{T,\ell}>20$ GeV, $-2.4 < \eta_\ell < 2.4$", fontsize=20)

pl.xlabel(r'$y$', fontsize=22, labelpad=-2)
pl.ylabel(r'$\mathcal{P}_{\rm PV}(Q,y,q_T)\; /\; \mathcal{P}(Q,y,q_T)\times 10^{6}$', fontsize=20, labelpad=-2)
pl.xlim(-2.4, 2.4)
pl.ticklabel_format(style="sci")
pl.ylim(-2, 2)
pl.xticks(size=20)
pl.yticks(size=20)

pl.plot(data[:10000,1], 1000000 * data[:10000,2],           'r-', label = '$q_T = 1$ GeV')
pl.plot(data[10000:20000,1], 1000000 * data[10000:20000,2], 'b-', label = '$q_T = 3$ GeV')
pl.plot(data[20000:30000,1], 1000000 * data[20000:30000,2], 'g-', label = '$q_T = 10$ GeV')

pl.legend(fontsize=20)

pl.savefig('PhaseSpaceRedFactor.pdf')
pl.show()

