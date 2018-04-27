#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

df=np.loadtxt("./ek_nphi_18_nel_6")
E0=df[:,1].min()

fig=plt.figure(figsize=(8,6))
ax=fig.add_subplot(111)
for i in range(6):
    ax.plot(df[:,0],df[:,i+1]-E0,marker='_',mfc='r',mec='r',ms=12,mew=2,ls='none')

ax.set_xlabel(r'$kl$',fontsize=22)
ax.set_ylabel(r'$\Delta E$',fontsize=22)
ax.set_xticks([0,0.5,1.0,1.5,2.0,2.5])
ax.set_xticklabels([0,0.5,1.0,1.5,2.0,2.5],fontsize=20)
ax.set_yticks([0.0,0.1])
ax.set_yticklabels([0.0,0.1],fontsize=20)
ax.set_xlim(-0.05,2.54)
ax.set_ylim(-0.001,0.194)
ax.text(0.2,0.03,r'$\nu=\frac{6}{18}$',fontsize=22)

plt.savefig("./excitation_spectrum.pdf")
plt.show()

