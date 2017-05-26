import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

beta = 100.0
U_list = [5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.3,6.5,6.7]
#U_list = [0.0]
symbol_list = ['-o','-^','-*','-s','-8','-h','-H','-s','-o','-o','-^','-s',"-o"]
for i,U in enumerate(U_list):
	d = loadtxt("Gf_tau_U%s_beta%s"%(U,beta))
	l = 0.0
        npoint = 1
	for i in range(npoint):
		l = l + d[-i-1,1]-d[-i-2,1]
        	#print l
	l = l/float(npoint)
	print U, d[-1,1]+ l/2, d[-1,1], npoint
		
