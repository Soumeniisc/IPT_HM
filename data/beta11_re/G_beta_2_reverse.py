import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t=1.0
beta=11.0
N_tau = 1024*2
dos = []
U_list = [1.0,2.0,3.0,4.0,5.0,6.0,7.0]
U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,5.0,5.5]
U_list = [0.0,0.4,1.0,1.8,2.2,3.0,3.4,4.2,4.6,5.0,5.5]

U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,4.7,4.8,4.9,5.0,5.2,5.1,5.3,5.4,5.5,5.6,5.7,5.9,6.1,6.3]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-o','-^','-s',"-o"]
data = open("dos_beta%s_re.dat"%beta,'w')
print >> data ,"# U/t, beta*G(beta/2)"
for i,U in enumerate(U_list):
	d = loadtxt("Gf_tau_U%s_beta%s_re"%(U,beta))
	dos.append(-1*d[N_tau/2][1])
        print >> data, U, beta*dos[i]
plt.plot(U_list,beta*np.array(dos),'-o',label = r'$\beta t =%s$'%beta)
plt.plot([U_list[0],U_list[-1]],[0,0],'-')
#plt.text(0.2,-0.5,r'$U_{mott}=6.0$',size=25)
plt.legend(loc=8)
plt.ylabel(r"$\beta Gf(\beta/2)$",size=20)
plt.xlabel("U/t",size=20)
plt.savefig("G_beta_2_HM_beta%s.eps"%beta,format="eps")
plt.show()
