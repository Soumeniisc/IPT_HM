import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

beta = 100.0
N=1025
e = -0.1
nf = 1.0/(1.0+exp(e*beta))


d1  = loadtxt("Gf_wn")

plt.plot(d1[:,0],d1[:,2],'-*',label='in_imag')
#plt.plot(wn,np.array(G).imag,'-*')
#plt.plot(wn,np.array(G).real,'-*',)
plt.legend()
plt.savefig("G_iw_imag.eps",format="eps")
plt.show()

d1  = loadtxt("Gf_tau")

plt.plot(d1[:,0],d1[:,1],'-*',label='in_imag')
#plt.plot(wn,np.array(G).imag,'-*')
#plt.plot(wn,np.array(G).real,'-*',)
plt.legend()
plt.savefig("G_tau.eps",format="eps")
plt.show()



