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

data = np.loadtxt('APFELvsTIMBA.dat')

pl.title(r"APFEL++ vs. TIMBA at $\mathcal{O}(\alpha_s)$", fontsize=22)

pl.xlabel('$p_T$ [GeV]', fontsize=22, labelpad=-2)
pl.ylabel('APFEL++ / TIMBA', fontsize=22, labelpad=-2)
pl.xlim(3.5, 12.5)
pl.ylim(0.98, 1.02)
pl.xticks(size=20)
pl.yticks(size=20)

ax.text(10, 1.015, r'$\frac{d\sigma}{dx dQ^2 dp_T d\eta}$', fontsize=33)

ax.text(4, 0.992, r'$\sqrt{S} = 300.9$ GeV', fontsize=22)
ax.text(4, 0.989, r'$x = 5\cdot 10^{-5}$', fontsize=22)
ax.text(4, 0.986, r'$Q^2 = 2$ GeV$^2$', fontsize=22)
ax.text(4, 0.983, r'$\eta = 1$', fontsize=22)

ax.text(8.5, 0.986,  r'PDF set: CT14nlo', fontsize=22)
ax.text(8.5, 0.983, r'FF set: DSS07 ($h+\overline{h}$)', fontsize=22)

pl.plot(data[:,2], data[:,5] / data[:,5], 'k-')
pl.plot(data[:,2], data[:,4] / data[:,5], 'r-')


pl.legend(fontsize=22)

pl.savefig('APFELvsTIMBA.pdf')
pl.show()
