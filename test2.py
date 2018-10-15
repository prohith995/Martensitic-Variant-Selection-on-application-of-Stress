## Date: 29/08/2016 / 02/10/2016
## Name: Rohith.P
## ??

import numpy as np
import math
from testdata2 import *


#def dotproduct(a,b):
#	dp = 0
#	for i in range(3):
#		dp = dp + a[i]*b[i]
#	return (dp)

def norm(a):
	nm = 0
	for i in range(len(a)):
		nm = nm + a[i]*a[i]
	return (math.sqrt(nm))

def unit_vector(a):
	uv = [x/norm(a) for x in a]
	return np.matrix(uv)

def costheta(a,b):
	ct = dotproduct(a,b) / (norm(a)*norm(b))
	return ct

def angle(a,b):
	angle = math.acos(costheta(a,b))
	return angle

#Q = np.matrix([357.6, 60.8, 4.7])
#sigma_gamma1 = np.dot(Q,stress)
#sigma_gamma = np.dot(sigma_gamma1 , Q.transpose())


# For a given stress state
sx = 5
sy = 0
sz = 0

stress = np.array([[sx,0,0],[0,sy,0],[0,0,sz]])

# For one particular bainite variant
p = [-0.168640, -0.760394, -0.627185]
d = [0.189920, -0.681934, 0.706326]
d1 = np.dot(d,p)
d2 = np.dot(d1,p)
e1 = d - d2 # d is resolved into delta*p and s*e where p and e are unit vectors
e = norm(e1) # making e as a unit vector

t = np.dot(stress, p)

sn = np.dot(t, p)
sigma_normal = np.dot(sn,p)

taumax = t - sigma_normal
tau = np.dot(taumax,e)

shear_stress = np.dot(tau, e)

