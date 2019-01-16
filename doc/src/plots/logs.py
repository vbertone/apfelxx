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

data = np.loadtxt('logs.dat')

pl.title(r"Non-perturbative modified logarithms at $Q = 10$ GeV", fontsize=20)

pl.xlabel(r'$q_T$ [GeV]', fontsize=22, labelpad=-2)
pl.ylabel(r'$I_p(q_T,Q)$', fontsize=22, labelpad=-2)
pl.xlim(0.01, 10)
pl.ylim(-2.5, 1)
pl.xticks(size=20)
pl.yticks(size=20)

#ax.text(10, 1.015, r'$\frac{d\sigma}{dx dQ^2 dp_T d\eta}$', fontsize=33)

#ax.text(4, 0.992, r'$\sqrt{S} = 300.9$ GeV', fontsize=22)
#ax.text(4, 0.989, r'$x = 5\cdot 10^{-5}$', fontsize=22)
#ax.text(4, 0.986, r'$Q^2 = 2$ GeV$^2$', fontsize=22)
#ax.text(4, 0.983, r'$\eta = 1$', fontsize=22)

#ax.text(8.5, 0.986,  r'PDF set: CT14nlo', fontsize=22)
#ax.text(8.5, 0.983, r'FF set: DSS07 ($h+\overline{h}$)', fontsize=22)

#pl.plot(data[:,2], data[:,5] / data[:,5], 'k-')
pl.plot(data[0:1000,1], data[0:1000,2], 'k-', label=r'$p=0$')
#pl.plot(data[0:1000,1], data[0:1000,4], 'k-.')

pl.plot(data[1000:2000,1], data[1000:2000,2], 'r-', label=r'$p=1$')
pl.plot(data[1000:2000,1], data[1000:2000,3], 'r--')
#pl.plot(data[1000:2000,1], data[1000:2000,4], 'r-.')

pl.plot(data[2000:3000,1], data[2000:3000,2], 'b-', label=r'$p=2$')
pl.plot(data[2000:3000,1], data[2000:3000,3], 'b--')
#pl.plot(data[2000:3000,1], data[2000:3000,4], 'b-.')

pl.plot(data[3000:4000,1], data[3000:4000,2], 'g-', label=r'$p=3$')
pl.plot(data[3000:4000,1], data[3000:4000,3], 'g--')
#pl.plot(data[3000:4000,1], data[3000:4000,4], 'g-.')

pl.legend(fontsize=22)

pl.savefig('logs.pdf')
pl.show()
