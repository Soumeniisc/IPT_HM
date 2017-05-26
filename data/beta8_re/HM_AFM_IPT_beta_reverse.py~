'''
to run this code u need following thing to make remeber
1> check the init para. if it is 1 then excutable read GB_old and GBd_old and density_old and densityd_old of previous U. other wise its read density_old and densityd_old  and  create GB and GBd of non interacting Green's function in HM on bethe lattice
2> if init=1 then turn on intialise(U,staggere). there U is the previous U. that GB_%U, GBd_%U, density_%U, densityd_%U would be copied to GB_old, GBd_old, density_old, densityd_old
3> chose U_i, for which U_i you want have calculation.
'''

from scipy import * 
import os,sys,subprocess
import numpy as np




def replace(old_file_,new_file_,stagger): #its put intial small staggerd magnetization for given U for A sublattice
	with open(new_file_, 'w') as new_file:
		with open(old_file_, 'r') as old_file:
			for i, line in enumerate(old_file):
				text = line.split()
				print text
				if i == 0: new_file.write("%d	%f	%f"%(int(text[0]),float(text[1])+stagger,float(text[2])-stagger))
				
def replaced(old_file_,new_file_,stagger): #its put intial small staggerd magnetization for given U for B sublattice
	with open(new_file_, 'w') as new_file:
		with open(old_file_, 'r') as old_file:
			for i, line in enumerate(old_file):
				text = line.split()
				if i == 0: new_file.write("%d	%f	%f"%(int(text[0]),float(text[1])-stagger,float(text[2])+stagger))
		#old_file.close()



def copy(input_file, copy_to):
	cmd = "cp %s  %s"%(input_file, copy_to)
	print os.popen(cmd).read()

def intialise(U,stagger):
	copy("density_U%s_beta%s_re"%(U,beta), "density_old")
	#copy("densityd_%s"%U, "densityd_old")
	#replace("density_%s"%U, "density_old",stagger)
	#replaced("densityd_%s"%U, "densityd_old",stagger)
	copy("Gf_U%s_beta%s_re"%(U,beta), "Gf_old")
	#copy("GfBd_%s"%U, "GfBd_old")

U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.9,6.1,6.3,6.5]
for i in range(len(U_list)-1):
        U = U_list[-2-i]
	# Parameters of the nonequilibrium DMFT simulation for the single-band Hubbard model
	dos="semicircular"  # density of states
	beta=10.0       # inverse temperature	
	U_i=U        # initial value of the interaction	
	N_tau=1024*2        # number of imaginary-time steps 2^11	
	N_iter=1000       # maximum iteration of DMFT self-consistency
	tolerance=0.0000001  # tolerance for DMFT convergence 10^-10
	solver="IPT"        # impurity solver
	mix=0.5		 # G = (1-mix)*G_new + G except first iteration
	delta=0.0	# orbital potential
        print "reading from:",U_list[-1-i],"  for U=", U_list[-2-i]
	if i==0: 
		init=   1        #0 when G0 is non-interacting GF.1 when G0 is interacting GF of previous U value
		stagger = 0.00
		intialise(U_list[-1-i],stagger)
	if i!=0: 
		init=1            #0 when G0 is non-interacting GF.1 when G0 is interacting GF of previous U value
		stagger = 0.00
		intialise(U_list[-1-i],stagger)

	cmd_ex = "./main_HM_eq.x %s %f  %f   %i %i %f %s %i %f %f"%(dos, beta,  U_i, N_tau, N_iter, tolerance, solver,init, mix, delta)
	print cmd_ex

	IPTA_info = open( "IPT_U%s_delta%s.dat"%(U_i,delta),'w')		
	subprocess.call(cmd_ex,shell=True,stdout=IPTA_info,stderr=IPTA_info)
	IPTA_info.flush()

	# coping for up spin
	copy("Gf", "Gf_U%s_beta%s_re"%(U_i,beta))
        copy("Gf_tau", "Gf_tau_U%s_beta%s_re"%(U_i,beta))
        #copy("G0", "G0_U%s_beta%s"%(U_i,beta))
        #copy("G0_tau", "G0_tau_U%s_beta%s"%(U_i,beta))
        #copy("Sig", "Sig_U%s_beta%s"%(U_i,beta))
        #copy("Sig_tau", "Sig_tau_U%s_beta%s"%(U_i,beta))
	copy("density", "density_U%s_beta%s_re"%(U_i,beta))
        copy("observables", "observables_U%s_beta%s_re"%(U_i,beta))
	
