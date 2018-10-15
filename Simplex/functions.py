
## This program contains all the useful and repetitive functions used in other files
import math
import numpy as np

# Calculating the unit vector of any vector
def unit_vector(a):
	uv = [x/norm(a) for x in a]
	return uv

# Calculating the norm of a vector
def norm(a):
	nm = 0
	for i in range(len(a)):
		nm = nm + a[i]*a[i]
	return (math.sqrt(nm))

# Calculated the shape deformation matrix for a given variant
def shape_deformation(m,d,p) :
	F = [ [1+m*d[0]*p[0], m*d[0]*p[1] , m*d[0]*p[2]],
		  [m*d[1]*p[0], 1+m*d[1]*p[1] , m*d[1]*p[2]],
		  [m*d[2]*p[0], m*d[2]*p[1] , 1+m*d[2]*p[2]]  ]
	return F

# Calculated the transformation matrix for a given euler angles
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

# Transforms any given tensor from smaple frame to crystal frame
def transform(Q,tensor):
	t1 = np.dot(Q,tensor)
	tf = np.dot(t1, Q.transpose())
	return tf

# Transforms any given tensor from crystal frame to sample frame
def transform_inverse(Q,tensor):
	t1 = np.dot(Q.transpose(),tensor)
	tf = np.dot(t1, Q)
	return tf

# Calculates the double dot product of 2 tensors
def double_dotproduct (a,b) :
	doubdot = 0
	for i in range(len(a)) :
		for j in range(len(b)) :
			doubdot = doubdot + a[i][j]*b[j][i] 
	return doubdot

# Converts a given strain into 5 dimensional deviatoric strain
def fiveD(i):
	j = [0]*5
	j[0] = (math.sqrt(3)+1)*i[1][1]*0.5 + (math.sqrt(3)-1)*i[2][2]*0.5
	j[1] = (math.sqrt(3)-1)*i[1][1]*0.5 + (math.sqrt(3)+1)*i[2][2]*0.5
	j[2] = math.sqrt(2)*i[1][2]
	j[3] = math.sqrt(2)*i[2][0]
	j[4] = math.sqrt(2)*i[0][1]
	return j

# Converts a given strain into 6 dimenbsional von mises strain
def sixD(i):
	j = [0]*6
	j[0] = i[0][0]
	j[1] = i[1][1]
	j[2] = i[2][2]
	j[3] = i[1][2]
	j[4] = i[0][2]
	j[5] = i[0][1]
	return j

# Obtain the 3x3 stress tensor from the 6 dimensional von mises strain
def sixD_inverse(i):
	s = [[0]*3 for j in range(3)]
	s[0][0] = i[0]
	s[1][1] = i[1]
	s[2][2] = i[2]
	s[1][2] = i[3]
	s[2][1] = i[3]
	s[0][2] = i[4]
	s[2][0] = i[4]
	s[0][1] = i[5]
	s[1][0] = i[5]
	return s

# Obtain the characteristic strain of a given variants using the crystallographic data using wang's method (newer method)
def strain_wang(m,d,p):
	F = np.asarray(shape_deformation( m, d, p))
	Ft = np.transpose(F)
	I = np.identity(3)
	dot = np.dot(Ft,F)
	strain = (dot - I)/2
	return strain

# Obtain the characteristic strain of a given variants using the crystallographic data using hilbert's method (older method)
def strain_hilbert(m,d,p):
	F = np.asarray(shape_deformation( m, d, p))
	Ft = np.transpose(F)
	I = np.identity(3)
	dot = (Ft+F)
	strain = (dot/2 - I)	
	return strain
