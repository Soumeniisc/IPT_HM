import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t=1.0

data = loadtxt("all_observable_50.0")
plt.plot(data[:,0],data[:,2],'-o',label = r'$\beta t =50.0$')
plt.legend(loc=8)
plt.ylabel(r"$d$",size=20)
plt.xlabel("U/t",size=20)
plt.savefig("double_occu_beta50.0.eps",format="eps")
plt.show()

