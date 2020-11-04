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
yv    = [0.1, 0.5, 0.7, 1.4, 2, 10]

plt.axvline(x = 1, ymin = 0, ymax = 10, linewidth = 1.5, color = 'k', linestyle='--')
plt.xlabel(r"$\kappa$")
plt.ylabel(r"$P^{{\rm NS},(0)}(\kappa,y)$")
plt.xlim(0.1, 10)
plt.ylim(-15, 15)
plt.xscale("log")

for y in yv:
    plt.plot(kappa[kappa > (y > 1) / y], [P0qq(k, y) for k in kappa[kappa > (y > 1) / y]], label = r"$y= " + str(y) + " $")
plt.legend()
plt.savefig('GPDP0.pdf')
plt.close()

######################################################

y      = np.linspace(0.1, 10, num = 500)
kappav = [0.1, 0.5, 0.9, 1.1, 2, 10]

plt.axhline(y = 0, xmin = 0, xmax = 10, linewidth = 1.5, color = 'k', linestyle='--')
plt.xlabel(r"$y$")
plt.ylabel(r"$P^{{\rm NS},(0)}(\kappa,y)$")
plt.xlim(0.1, 10)
plt.ylim(-20, 10)
plt.xscale("log")

for kappa in kappav:
    plt.plot(y[y < 10 - 9 * (kappa < 1)], [P0qq(kappa, ys) for ys in y[y < 10 - 9 * (kappa < 1)]], label = r"$\kappa= " + str(kappa) + " $")
plt.legend()
plt.savefig('GPDP02.pdf')
plt.close()




