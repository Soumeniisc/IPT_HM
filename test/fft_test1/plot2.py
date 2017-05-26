import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

beta = 100.0
N=1025
e = -0.1

tau = []
N_t = 100
G_t = []
for i in range(N_t):
        tau.append(i*beta/N_t)
	G_t.append(-exp(-e*tau[i])/(1+exp(-e*beta)))
d  = loadtxt("G_tau_free")
#d1 = loadtxt("G_tau.dat")


#plt.plot(tau,G_t,label='-*')
plt.plot(d[:,0],d[:,1],'-o',label='fortran')
plt.plot(tau,G_t,'-*')
#plt.plot(d1[:,0],d1[:,1],'-*',label='python')
plt.legend()
plt.show()

wn=[]
G=[]
for i in range(N):
	wn.append((2*i+1.0)*np.pi/beta)	
	G.append(1.0/(1j*wn[i] -e))

d  = loadtxt("G_wn_out")
plt.plot(d[:,0],d[:,1],'-o',label='real')
plt.plot(d[:,0],d[:,2],'-o',label='imag')
plt.plot(wn,np.array(G).imag,'-*')
plt.plot(wn,np.array(G).real,'-*',)
plt.legend()
plt.show()

d  = loadtxt("G_wn_out")
d1  = loadtxt("G_wn")
plt.plot(d[:,0],d[:,2],'-o',label='out_imag')
plt.plot(d1[:,0],d1[:,2],'-*',label='in_imag')
#plt.plot(wn,np.array(G).imag,'-*')
#plt.plot(wn,np.array(G).real,'-*',)
plt.legend()
plt.show()


