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

data = np.loadtxt('LeptonicModulation.dat')

pl.title(r"Leptonic modulation at $Q = M_Z$ and $y = 0$", fontsize=20)

pl.xlabel(r'$p_{T,\ell}$ [GeV]', fontsize=22, labelpad=-2)
pl.ylabel(r'$\frac{d\sigma}{dQ dy dq_T dp_T} / \frac{d\sigma}{dQ dy dq_T}\;\times 100$', fontsize=22, labelpad=-2)
pl.xlim(41.5, 45.5)
pl.ylim(0, 20)
pl.xticks(size=20)
pl.yticks(size=20)

np = 10000

pl.plot(data[0*np:1*np,1], 100*data[0*np:1*np,2], 'r-', label=r'$q_T = 1$ GeV')
pl.plot(data[1*np:2*np,1], 100*data[1*np:2*np,2], 'b-', label=r'$q_T = 2$ GeV')
pl.plot(data[2*np:3*np,1], 100*data[2*np:3*np,2], 'm-', label=r'$q_T = 3$ GeV')
pl.plot(data[3*np:4*np,1], 100*data[3*np:4*np,2], 'g-', label=r'$q_T = 4$ GeV')
pl.plot(data[4*np:5*np,1], 100*data[4*np:5*np,2], 'c-', label=r'$q_T = 5$ GeV')

pl.legend(fontsize=22)

pl.savefig('LeptonicModulation.pdf')
pl.show()
