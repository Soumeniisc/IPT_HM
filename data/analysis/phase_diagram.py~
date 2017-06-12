import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

beta_list = [100.0, 50.0, 25.0, 20.0, 15.0, 12.0, 11.0]
U_increase =[6.7,   6.5,  6.0,  5.8,  5.4,  5.1,  5.0]
U_decrease =[5.2,   5.2,  5.2,  5.2,  5.1,  5.0,  5.0]


plt.plot(U_increase,1.0/np.array(beta_list),'-o', label="U_increase")
plt.plot(U_decrease,1.0/np.array(beta_list),'-o', label="U_decrease")
plt.legend()
plt.ylabel(r"$T/t$",size=20)
plt.xlabel("U/t",size=20)
plt.xlim(4,7)
plt.ylim(0.0,0.1)
plt.savefig("HM_paramagnetic_phase_diagram.eps",format="eps")
plt.show()
