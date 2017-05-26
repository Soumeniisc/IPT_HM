import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t=1.0

beta_list = [10.0,11.0,12.0,15.0,20.0,25.0,50.0,100.0]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-o','-^','-s',"-o"]
for i,beta in enumerate(beta_list):
        data = loadtxt("dos_beta%s.dat"%beta)
	plt.plot(data[:,0],data[:,1],symbol_list[i],label = r'$\beta t =%s$'%beta)
plt.legend(loc=8)
plt.ylabel(r"$\beta Gf(\beta/2)$",size=20)
plt.xlabel("U/t",size=20)
plt.savefig("dos_0_HM_diff_beta.eps",format="eps")
plt.show()
