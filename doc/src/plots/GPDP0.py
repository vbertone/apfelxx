import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

def P0qq(kappa, y):
    CF = 4. / 3.
    if kappa < 1:
        return 2 * CF * ( 1 + ( 1 - 2 * kappa**2 ) * y**2 ) / ( ( 1 - y ) * ( 1 - kappa**2 * y**2 ) )
    else:
        return 2 * CF * ( 1 / ( 1 - y ) + ( 1 - kappa ) / 2 / ( kappa * ( 1 + kappa * y ) ) )

kappa = np.linspace(0, 10, num = 500)
yv    = [0.01, 0.1, 0.3, 0.5, 0.7]
P0    = [[P0qq(k, y) for k in kappa] for y in yv]

plt.axvline(x = 1, ymin = 0, ymax = 10, linewidth = 1.5, color = 'k', linestyle='--')
plt.xlabel(r"$\kappa$")
plt.ylabel(r"$P^{{\rm NS},(0)}(\kappa,y)$")
plt.xlim(0.1, 10)
plt.ylim(1, 16)
plt.xscale("log")
#plt.yscale("log")

for i in range(len(P0)):
    plt.plot(kappa, P0[i], label = r"$y= " + str(yv[i]) + " $")
plt.legend()

plt.savefig('GPDP0.pdf')
