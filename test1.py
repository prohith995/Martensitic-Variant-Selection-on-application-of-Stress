## Date: 23/08/2016
## Name: Rohith.P
## BTP first task

import random
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

def costheta(a,b):
	ct = dotproduct(a,b) / (norm(a)*norm(b))
	return ct

a = [1,2,3]
b = [2,-1,3]
#print costheta(a,b)

stress = [2,5,1]
stress_norm = norm(stress)
stress_dir = [stress[i]/norm(stress) for i in range(3)]

#print stress_norm , stress_dir

m = stress_norm

def interaction_energy(hp, disp):

	global stress_norm
	global stress
	global m

	sigma_n = stress_norm * costheta(stress, hp)
	tau = stress_norm * math.sqrt(1-costheta(stress, hp)**2)
	delta = m * dotproduct(disp, hp)
	s = math.sqrt(m**2 - delta**2)

	eng = sigma_n*delta + tau*s
	
	return(eng)

# Calculating the interaction energy for the different variants

inteng = [[0]*len(hpnormal) for i in range(2)]
inteng[1] = [1,2,3,4,5,6]

for n in range(len(hpnormal)):
	print hpnormal[n], dispdir[n]
	inteng[0][n] = interaction_energy( hpnormal[n], dispdir[n])
	

s1 = random.uniform(0,1)
s2 = random.uniform(0,1)
s3 = random.uniform(0,1)


# Insertion Sort
for i in range(0, len(hpnormal)):
    temp = inteng[0][i]
    temp1 = inteng[1][i]
    j = i - 1
    while j >= 0 and inteng[0][j] > temp : 
        inteng[0][j+1] = inteng[0][j]
	inteng[1][j+1] = inteng[1][j]
        j = j - 1
    
    inteng[0][j+1] = temp
    inteng[1][j+1] = temp1

print inteng