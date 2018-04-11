#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

df_nel4=np.loadtxt("./gs_energy_nel4")
df_nel5=np.loadtxt("./gs_energy_nel5")
df_nel6=np.loadtxt("./gs_energy_nel6")

fig=plt.figure(figsize=(8,8))
ax=fig.add_subplot(111)

ax.plot(df_nel4[:,0],df_nel4[:,1],ls='-',marker='o',mfc='r',mec='r',c='r',ms=8,label=r'$n=4$')
ax.plot(df_nel5[:,0],df_nel5[:,1],ls='none',marker='s',mfc='none',mec='b',mew=1,c='b',ms=8,label=r'$n=5$')
ax.plot(df_nel6[:,0],df_nel6[:,1],ls='none',marker='o',mfc='none',mec='k',mew=1,c='k',ms=8,label=r'$n=6$')

ax.set_xlim(0.245,0.6)
ax.set_yticks([-0.5,-0.45,-0.4,-0.35])
ax.set_yticklabels([-0.50,-0.45,-0.40,-0.35],fontsize=14)
ax.set_xticks([0.30,0.40,0.50,0.60])
ax.set_xticklabels([0.30,0.40,0.50,0.60],fontsize=14)
ax.set_ylim(-0.5,-0.348)
ax.set_xlabel(r'$\nu$',fontsize=22)
ax.set_ylabel(r'$E_{gs}/(e^2/\epsilon l)$',fontsize=22)
ax.axvline(x=0.33333,ls='--',c='k')
ax.axvline(x=0.4,ls='--',c='k')
ax.axvline(x=0.5,ls='--',c='k')

plt.legend(loc='best',numpoints=1,fontsize=20)
plt.savefig("./ground_state_energy_vs_nu.pdf")
plt.show()
