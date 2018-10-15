
## Penrose inverse

import numpy as np
import math
import functions as fn
import data_variants as data
from operator import itemgetter

# This function is used to calculate the poles from the given preferred variants
def obtain_poles(pref_var,euler_angles):
	# Obtaining the normal directions of all the variants in the preferred variants list
	if pref_var == None :
		return []
	norm_pref_var1 = [data.p[int(i)-1] for i in pref_var] 	
	norm_pref_var = []
	for i in norm_pref_var1:
		norm_pref_var.append(i)
		norm_pref_var.append([-i[0],-i[1],-i[2]])
	
	# Obtaining the normal vectors in the sample frame
	Q = np.asarray(fn.transformation(euler_angles))
	norm_pref_var_samp_fram = [fn.transform_inverse(Q,i) for i in norm_pref_var]
	
	poles = []
	for i in norm_pref_var_samp_fram :
		if i[2] >= 0 : ## To see if the pole exists in the pole figure. i.e. it is in the northen hemisphere
			x = i[0]/(i[2]+1)
			y = i[1]/(i[2]+1)
			poles.append([x,y])
	
	return poles

# This function solves the normal equation to give the values of the coefficients
def solve_normal(strain) :

	Mt = np.asarray(fn.M(6))
	M = np.transpose(Mt)
	m1 = np.dot(Mt,M)
	inv = np.linalg.pinv(m1, rcond=1e-15)
	m2 = np.dot(inv,Mt)

	X = np.dot(m2,strain)
	return X

# This fucntion solves the normal equation for all the grains and gives the preferred variants and other required outputs
def solve(n,strain):

	poles = []
	ener = []
	c = 0
	prefvars = []

	f = open("rand5000.txt",'r')
	#f = open("fcc_roll.txt",'r')
	for i in range(n):
		variants = []
		prefvar = []
		a = f.readline().split()
		euler_angle = [float(a[0]), float(a[1]), float(a[2])]
		print "Euler Angle No. " , i+1, ": ", euler_angle
		
		Q = np.asarray(fn.transformation(euler_angle))
		s1 = fn.sixD(fn.transform(Q,strain))
		s = np.transpose(s1)
		x = solve_normal(s)	

		for i in range(len(x)) :
			if x[i]>0: 
				variants.append([i+1, x[i]])

		variants = sorted(variants, key=itemgetter(1), reverse = True)
		
		maxener = variants[0][1]
		for i in variants :
			if i[1] >= 0.1*maxener :
				prefvar.append(i[0])
				ener.append(i[1])
		prefvars.extend(prefvar)

		poles.extend(obtain_poles(prefvar,euler_angle))

	return (poles, prefvars, ener)
			