import matplotlib.pyplot as plt

#latex font and text-rendering
plt.rc('font', family='serif')
plt.rc('mathtext', **{'default': 'regular'})
plt.rc('text', usetex=True)

#figure style
plt.style.use('seaborn-v0_8-paper')
subplots_adjust = {'left':0.05, 'right':0.95, 'top':0.95, 'bottom':0.05,'wspace':0.20,'hspace':0.40}
plt.rc('xtick', **{'direction': 'in'})
plt.rc('ytick', **{'direction': 'in'})
