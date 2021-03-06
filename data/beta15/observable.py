import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t=1.0
beta=15.0
U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,4.6,4.7,4.8,4.9,5.0,5.2,5.3,5.4,5.5,5.6,5.7,5.9,6.1,6.3,6.5]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-o','-^','-s',"-o"]
data = open("all_observable_%s"%beta,'w')
print >> data, " # U_i(1)   density(2)   double_occu(3)  kinetic_energy(4) interaction_enegy(5) total_energy(6)"
for i,U in enumerate(U_list):
	d = loadtxt("observables_U%s_beta%s"%(U,beta))
	#print d
	print >> data, d[0], d[1], d[2], d[3], d[4], d[5]
data.close()

data = loadtxt("all_observable_%s"%beta)
plt.plot(data[:,0],np.array(data[:,2]),"-o")
plt.text(3.5,0.1,r'$\beta=%s$'%beta,size=20)
plt.legend()
plt.ylabel("double occu")
plt.xlabel("U/t")
plt.savefig("double_occu_U_HM_beta%s.eps"%beta,format="eps")
plt.show()

plt.plot(data[:,0],np.array(data[:,1]),"-o")
plt.text(3.5,0.5,r'$\beta=%s$'%beta,size=20)
plt.legend()
plt.ylabel("density")
plt.xlabel("U/t")
plt.savefig("density_U_HM_beta%s_re.eps"%beta,format="eps")
plt.show()

plt.plot(data[:,0],np.array(data[:,3]),"-o", label="kinetic")
plt.plot(data[:,0],np.array(data[:,4]),"-o", label="potential")
plt.plot(data[:,0],np.array(data[:,5]),"-o", label="total")
plt.text(3.5,0.0,r'$\beta=%s$'%beta,size=20)
plt.legend()
plt.ylabel("energy")
plt.xlabel("U/t")
plt.savefig("double_occu_U_HM_beta%s.eps"%beta,format="eps")
plt.show()

