## Date: 29/08/2016 / 02/10/2016
## Name: Rohith.P
## Veryfying for 24 variants : Pg 57 of full thesis
# Interaction energy for all 24 variants calculated using wangs method and hilbert's method

import matplotlib.pyplot as plt
import pylab
import numpy as np
import math
from test1data import *

def dotproduct(a,b):
	dp = 0
	for i in range(3):
		dp = dp + a[i]*b[i]
	return (dp)

def norm(a):
	nm = 0
	for i in range(len(a)):
		nm = nm + a[i]*a[i]
	return (math.sqrt(nm))

def unit_vector(a):
	uv = [x/norm(a) for x in a]
	return uv

def costheta(a,b):
	ct = dotproduct(a,b) / (norm(a)*norm(b))
	return ct

def angle(a,b):
	angle = math.acos(costheta(a,b))
	return angle

def transformation(av):
	(a1,a,a2) = (z for z in av)
	ss = [[0]*3 for x in range(3)]
	ss[0][0] = math.cos(math.radians(a1))*math.cos(math.radians(a2))-math.sin(math.radians(a1))*math.cos(math.radians(a))*math.sin(math.radians(a2))
	ss[0][1] = math.sin(math.radians(a1))*math.cos(math.radians(a2))+math.cos(math.radians(a1))*math.cos(math.radians(a))*math.sin(math.radians(a2))
	ss[0][2] = math.sin(math.radians(a))*math.sin(math.radians(a2))
	ss[1][0] = -math.cos(math.radians(a1))*math.sin(math.radians(a2))-math.sin(math.radians(a1))*math.cos(math.radians(a))*math.cos(math.radians(a2))
	ss[1][1] = -math.sin(math.radians(a1))*math.sin(math.radians(a2))+math.cos(math.radians(a1))*math.cos(math.radians(a))*math.cos(math.radians(a2))
	ss[1][2] = math.sin(math.radians(a))*math.cos(math.radians(a2))
	ss[2][0] = math.sin(math.radians(a1))*math.sin(math.radians(a))
	ss[2][1] = -math.cos(math.radians(a1))*math.sin(math.radians(a))
	ss[2][2] = math.cos(math.radians(a))
	return ss


def interaction_energy (stress, strain) :
	doubdot = 0
	for i in range(len(stress)) :
		for j in range(len(strain)) :
			doubdot = doubdot + stress[i][j]*strain[j][i] 

	return doubdot

def shape_deformation(m,d,p) :
	F = [ [1+m*d[0]*p[0], m*d[0]*p[1] , m*d[0]*p[2]],
		  [m*d[1]*p[0], 1+m*d[1]*p[1] , m*d[1]*p[2]],
		  [m*d[2]*p[0], m*d[2]*p[1] , 1+m*d[2]*p[2]]  ]
	return F

def strain_wang(m,d,p):
	F = np.asarray(shape_deformation( m, d, p))
	Ft = np.transpose(F)
	I = np.identity(3)
	dot = np.dot(Ft,F)
	strain = (dot - I)/2
	return strain

def strain_hilbert(m,d,p):
	F = np.asarray(shape_deformation( m, d, p))
	Ft = np.transpose(F)
	I = np.identity(3)
	dot = (Ft+F)
	strain = (dot/2 - I)	
	return strain

def transform(Q,tensor):
	t1 = np.dot(Q,tensor)
	tf = np.dot(t1, Q.transpose())
	return tf



def variant_selection_wang():

	for i in range(len(p1)):

		strain_wan = strain_wang(m,unit_vector(d1[i]), unit_vector(p1[i]))
		U[i] = interaction_energy(stress,strain_wan)
	
	U1 = [0]*len(U)

	for i in range(len(U1)):
		U1[i] = [i+1, U[i]/math.pow(10,6)]
	return U1


def variant_selection_hilbert():

	V = [0]*len(p1)
	for i in range(len(p1)):
		strain_hil = strain_hilbert(m,unit_vector(d1[i]), unit_vector(p1[i]))
		V[i] = interaction_energy(stress,strain_hil)
	V1 = [0]*len(V)
	for i in range(len(V1)):
		V1[i] = [i+1, V[i]/math.pow(10,6)]

	return V1


def sorting (A):
	for i in range(len(A)) :
		for j in range(1,len(A)-i):
			if A[j][1] > A[j-1][1] :
				temp = A[j]
				A[j] = A[j-1]
				A[j-1] = temp
	return A


def macroscopic_strain_hilbert(n,v):

	Q = np.asarray(transformation(euler_angles))
	trans_strain = np.asarray([[0,0,0],[0,0,0],[0,0,0]])
	weight = 0.63/n
	for i in range(n):
		index = v[i][0]-1
		strain = strain_hilbert(m,unit_vector(d1[index]),unit_vector(p1[index]))
		s_crystal = transform(Q.transpose(), strain)

		trans_strain = trans_strain + weight*s_crystal

	return trans_strain

def macroscopic_strain_wang(n,u):

	Q = np.asarray(transformation(euler_angles))
	trans_strain = np.asarray([[0,0,0],[0,0,0],[0,0,0]])
	weight = 0.63/n
	for i in range(n):
		index = u[i][0]-1
		strain = strain_wang(m,unit_vector(d1[index]),unit_vector(p1[index]))
		s_crystal = transform(Q.transpose(), strain)

		trans_strain = trans_strain + weight*s_crystal

	return trans_strain


def stress_calc():
	wang_ms = [0]*24
	hilbert_ms = [0]*24
	for i in range(1,25):
		wang_ms[i-1] = macroscopic_strain_wang(i,U)
		hilbert_ms[i-1] = macroscopic_strain_hilbert(i,V)

	return (wang_ms, hilbert_ms)

def plotting(wang_ms, hilbert_ms):

	longitudinal_strain_wang = [0]*24
	lateral_strain_wang1 = [0]*24
	lateral_strain_wang2 = [0]*24

	for i in range(24):
		longitudinal_strain_wang[i] = wang_ms[i][0][0]
		lateral_strain_wang1[i] = wang_ms[i][1][1]
		lateral_strain_wang2[i] = wang_ms[i][2][2]


	longitudinal_strain_hilbert = [0]*24
	lateral_strain_hilbert1 = [0]*24
	lateral_strain_hilbert2 = [0]*24

	for i in range(24):
		longitudinal_strain_hilbert[i] = hilbert_ms[i][0][0]
		lateral_strain_hilbert1[i] = hilbert_ms[i][1][1]
		lateral_strain_hilbert2[i] = hilbert_ms[i][2][2]

	active_variants = [(i+1) for i in range(24)]


	plt.plot(active_variants, longitudinal_strain_wang, label = "long_wang")
	plt.plot(active_variants, lateral_strain_wang1, label = "long_wang")
	plt.plot(active_variants, lateral_strain_wang2, label = "long_wang")	

	plt.plot(active_variants, longitudinal_strain_hilbert, label = "long_hil")
	plt.plot(active_variants, lateral_strain_hilbert1, label = "lat1_hil")
	plt.plot(active_variants, lateral_strain_hilbert2, label = "lat2_hil")	

	plt.legend()
	plt.title('Strain vs Active Variants')
	plt.xlabel('no of active variants')
	plt.ylabel('true strain')

	plt.show()



def run():
	
	global stress
	global U
	global V
	global Q
	global wang_ms
	global hilbert_ms

	w = 0
	h = 0

	Q = np.asarray(transformation(euler_angles))
	stress = transform(Q,initial_stress())		
	U = [0]*len(p1)
	V = [0]*len(p1)

	U = variant_selection_wang()
	V = variant_selection_hilbert()

	U1 = sorting (U)
	V1 = sorting (V)

	for i in range(len(U)):
		#print U1[i][0] , "  ", U1[i][1], "     ", V1[i][0],"  ", V1[i][1]

	for i in U1:
		if i[1] > 0 : w = w + 1
	for j in V1:
		if j[1] > 0 : h = h + 1
		
	return (w,h)



euler_angles = [128.99, 48.27, 24.27]
m = 0.229367


# For different types of stresses

def initial_stress():
	sx = 200*math.pow(10,6)
	stress_initial = np.array([[sx,0,0],[0,0,0],[0,0,0]])
#	stress_initial = np.array([[sx,0,0],[0,sx,0],[0,0,0]])
#	stress_initial = np.array([[sx,0,0],[0,sx,0],[0,0,sx]])
#	stress_initial = np.array([[-sx,0,0],[0,0,0],[0,0,0]])
#	stress_initial = np.array([[-sx,0,0],[0,-sx,0],[0,0,0]])
#	stress_initial = np.array([[-sx,0,0],[0,-sx,0],[0,0,-sx]])
#	stress_initial = np.array([[0,sx,0],[sx,0,0],[0,0,0]])
#	stress_initial = np.array([[0,-sx,0],[-sx,0,0],[0,0,0]])
#	stress_initial = np.array([[sx,0,0],[0,-sx,0],[0,0,-sx]])
#	stress_initial = np.array([[sx,0,0],[0,sx,0],[0,0,-sx]])

	return stress_initial

stress = initial_stress()
(w,h) = run()
wang1 = w
hilbert1 = h

#print wang1, hilbert1
#plotting(wang, hilbert)
