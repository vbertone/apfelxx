import numpy as np
import pandas as pd
import tests.ACOT_tests.matplotlib_settings as plt_settings
import matplotlib.pyplot as plt

path = "/home/peter/work/codes/apfelxx_ACOT/tests/ACOT_tests"
file_nCTEQ = path + "/nCTEQ_results/NC_SACOT-chi_test.csv"
file_apfel = path + "/APFELxx_results/NC_F2_SACOT-chi.csv"

def ij_to_csv_label(i,j):
  cur = ["NC","WP","WM"]
  strucf = ["F2","FL","F3"]

  return cur[i]+strucf[j]

colors = ["orange","red","royalblue","darkgreen"]

currents = [r"$\gamma/Z$"]
sfs = [r"$F_2$",r"$F_L$",r"$F_3$"]

#prepare data
nCTEQ = pd.read_csv(file_nCTEQ)
apfel = pd.read_csv(file_apfel)

Q2 = nCTEQ[nCTEQ["x"]==nCTEQ["x"].iloc[0]]["Q2"].to_numpy()
x = nCTEQ[apfel["Q2"]==Q2[0]]["x"].to_numpy()

#plot each structure function
fig,ax = plt.subplots(6,3,figsize=(12,9*6/4),gridspec_kw={'height_ratios': [2, 1, 2, 1, 2, 1]})
for i,current in enumerate(currents):
  for j,sf in enumerate(sfs):
    tot = ax[2*i][j]
    ratio = ax[2*i+1][j]
    label = ij_to_csv_label(i,j)
    for k,Q2k in enumerate(Q2):
      # if Q2i <= 1.5**2: continue
      
      val_nCTEQ = nCTEQ[nCTEQ["Q2"]==Q2k][label].to_numpy()
      val_apfel = apfel[apfel["Q2"]==Q2k][label].to_numpy()
      tot.plot(x,val_nCTEQ,color=colors[k],label='Q={Q:.2f}'.format(Q=np.sqrt(Q2k)))
      tot.plot(x,val_apfel,ls='dashed',color=colors[k])
      ratio.plot(x,val_apfel/val_nCTEQ-1,color=colors[k])
      # print("Q:",np.sqrt(Q2i))
      # print(val_apfel/val_nCTEQ)

    tot.legend()
    tot.set_yscale('log')
    #symmetric ylim around 1
    max_deviation = np.nanmax(abs(np.ma.masked_invalid(apfel[label]/nCTEQ[label])-1))
    print(label+':','max_deviation:',max_deviation)
    #clip plotting range if deviation is below 10^-4
    yabs_max = max(max_deviation*1.2,1e-4) 
    ratio.set_ylim(ymin=-yabs_max, ymax=yabs_max)
    ratio.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False)
    # ratio.ticklabel_format(useOffset=False)


    for plot in [tot,ratio]:
      plot.set_xscale('log')
      plot.set_xlim(min(x),max(x))

    tot.set_title(current+' '+sf)
    ratio.set_title(r"$\texttt{apfelxx}/\texttt{nCTEQ}-1$")
    # plt.suptitle(key)
plt.subplots_adjust(**plt_settings.subplots_adjust)
plt.show()  

#collect as PDF
#TODO