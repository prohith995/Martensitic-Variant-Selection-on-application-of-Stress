
import data_variants as data
import functions as fn
import numpy as np
import math
from operator import itemgetter

# This function calculates the characteristic strains of each variants
def M(n):
	M1 = []
	if n == 6 :
		for i in range(len(data.p)):
			s = fn.strain_wang( data.m, fn.unit_vector(data.d[i]), fn.unit_vector(data.p[i]))
			M1.append( fn.sixD(s) )
	else :
		for i in range(len(data.p)):
			s = fn.strain_wang( data.m, fn.unit_vector(data.d[i]), fn.unit_vector(data.p[i]))
			M1.append( fn.fiveD(s) )
	return M1

# Calculates the objective function required for the simplex table
def simplexing (euler_angle, strain) :
	Q = np.asarray(fn.transformation(euler_angle))
	s = fn.sixD(fn.transform(Q,strain))

	obj_fn = []
	for j in s:
		obj_fn.append(-1*j)
	return obj_fn

# This fucntion calls the simplex algorthm by giving the required inputs
def optimize_simplex (euler_angle, strain):
	
	import simplex
	constraints = []
	for i in M(6) :
		row = []
		for j in i :
			row.append(j)
			#row.append(-1*j)
		constraints.append([row,1])

	objective_function = simplexing(euler_angle, strain)
	table1 = simplex.table(objective_function)
	
	(variants, energy) = table1.run(constraints)
	return (variants, energy)

#	table = table1.run(constraints)
#	return table

# This is to obtain the poles from the given active variants
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

# This is the main function which runs the simplex algorithm for all the grains
def simplex_run(n,strain):

	poles = []
	energy = []
	prefvars = []

	f = open("rand5000.txt",'r')
	#f = open("fcc_roll.txt",'r')

	for i in range(n):
		a = f.readline().split()
		euler_angle = [float(a[0]), float(a[1]), float(a[2])]
		print "\n\nEuler Angle No. " , i+1, ": ", euler_angle
		
		(pref_var, ener) = optimize_simplex(euler_angle,strain)
		energy.append(ener)
		prefvars.extend(pref_var)
		poles.extend(obtain_poles(pref_var,euler_angle))

	return (poles, prefvars, energy)

# This obtains the final simplex table for the verification of stress
def simplex_obtain_stress(n,strain):

	poles = []
	variants = []
	energy = []

	f = open("rand5000.txt",'r')
	for i in range(n):
		a = f.readline().split()
		euler_angle = [float(a[0]), float(a[1]), float(a[2])]
	print "Euler Angle No. " , i+1, ": ", euler_angle
		
	table = optimize_simplex(euler_angle,strain)		
	return table

	