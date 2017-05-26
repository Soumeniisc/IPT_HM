import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t=1.0
beta=15.0
U_list = [1.0,2.0,3.0,4.0,5.0,6.0,7.0]
U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,5.0,5.5]
#U_list = [0.0,0.4,1.0,1.8]
U_list = [4.6,5.0,5.2,5.3,5.4,5.5,5.6,5.7]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-o','-^','-s',"-o"]
for i,U in enumerate(U_list):
	d = loadtxt("Gf_tau_U%s_beta%s"%(U,beta))
	plt.plot(d[:,0],np.array(d[:,1]),symbol_list[i],label = r'$U=%s$'%U)
	plt.text(0.2,-0.5,r'$U_{mott}=6.0$',size=25)
	plt.legend()
	plt.ylabel("Gf_tau")
	plt.xlabel("tau")
	#plt.xlim(0,10)
	#plt.ylim(-1 ,0)
plt.savefig("Gf_tau_HM_beta%s.eps"%beta,format="eps")
plt.show()

for i,U in enumerate(U_list):
	d = loadtxt("Gf_tau_U%s_beta%s"%(U,beta))
	plt.plot(d[:,0],np.array(d[:,1]),symbol_list[i],label = r'$U=%s$'%U)
	plt.text(0.2,-0.5,r'$U_{mott}=6.0$',size=25)
	plt.legend()
	plt.ylabel("Gf_tau")
	plt.xlabel("tau")
	plt.xlim(0,0.5)
	#plt.ylim(-1 ,0)
plt.savefig("Gf_tau_zoom_left_HM_beta%s.eps"%beta,format="eps")
plt.show()

for i,U in enumerate(U_list):
	d = loadtxt("Gf_tau_U%s_beta%s"%(U,beta))
	plt.plot(d[:,0],np.array(d[:,1]),symbol_list[i],label = r'$U=%s$'%U)
	plt.text(0.2,-0.5,r'$U_{mott}=6.0$',size=25)
	plt.legend()
	plt.ylabel("Gf_tau")
	plt.xlabel("wn")
	plt.xlim(beta-0.5,beta)
	#plt.ylim(-1 ,0)
plt.savefig("Gf_tau_zoom_right_HM_beta%s.eps"%beta,format="eps")
plt.show()

t=1.0


#U_list = [0.0,1.8,3.4,5.0,5.5,5.9,6.0,6.1,6.3,6.5]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-^','-s',"-o"]
for i,U in enumerate(U_list):
	d = loadtxt("Gf_U%s_beta%s"%(U,beta))
        # G0^{-1} = iwn-G(iwn)
	#G0= 1.0/( 1j*np.array(d[:,0]) -np.array(d[:,1]) -1j*np.array(d[:,2]) )
	plt.plot(d[:,0],np.array(d[:,2]),symbol_list[i],label = r'$U=%s$'%U)
	#plt.text(0.0,1.0,r'$Ui=0.0$',size=25)
	
	plt.ylabel("ImGf")
	plt.xlabel("wn")
	plt.xlim(0,5)
	plt.ylim(-1 ,0)
#plt.plot(d[:,0],-1.0/np.array(d[:,0]),'-s',label = r'$1/iwn$')
plt.legend()
plt.savefig("IMG_Gf_HM_beta%s.eps"%beta,format="eps")
plt.show()

for i,U in enumerate(U_list):
	d = loadtxt("Gf_U%s_beta%s"%(U,beta))
        # G0^{-1} = iwn-G(iwn)
	#G0= 1.0/( 1j*np.array(d[:,0]) -np.array(d[:,1]) -1j*np.array(d[:,2]) )
	plt.plot(d[:,0],np.array(d[:,1]),symbol_list[i],label = r'$U=%s$'%U)
	#plt.text(0.0,1.0,r'$Ui=0.0$',size=25)
	
	plt.ylabel("ReGf")
	plt.xlabel("wn")
	plt.xlim(0,5)
	#plt.ylim(-1 ,0)
#plt.plot(d[:,0],-1.0/np.array(d[:,0]),'-s',label = r'$1/iwn$')
plt.legend()
plt.savefig("Re_Gf_HM_beta%s.eps"%beta,format="eps")
plt.show()
