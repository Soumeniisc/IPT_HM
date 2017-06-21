import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


U_list = [1.0,3.0,5.0,6.5]
beta = 100.0
for U in U_list:
	data = loadtxt("DOS_U%s_beta100.0"%U)
	plt.plot(data[:,0],data[:,1],'-',label="U=%st"%U)

plt.legend()
plt.xlabel(r"$w$",size=20)
plt.xlabel(r"$A(w)$",size=20)
plt.text(-2.0,0.3,r"$\beta*t=%s$"%beta, size=20)
plt.savefig("dos_M_MI_transition_HM.eps", format="eps")
plt.show()


