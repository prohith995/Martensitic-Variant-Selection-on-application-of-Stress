## Date: 29/08/2016
## Name: Rohith.P
## BTP first task

import numpy as np
import math
from test1data import *

def norm(a):
	nm = 0
	for i in range(len(a)):
		nm = nm + a[i]*a[i]
	return (math.sqrt(nm))

def interaction_energy(p,d):
	global m
	global s
	u = m*norm(d)*norm(p)*s
	return u

m = 0.009981
s = 100

inteng = [[0]*(len(hpnormal)+1) for i in range(2)]
#print len(hpnormal)
for i in range(1,len(hpnormal)+1):
	inteng[1][i] = i

for n in range(len(hpnormal)):
#
	print hpnormal[n], dispdir[n]
	inteng[0][n] = interaction_energy( hpnormal[n], dispdir[n])
	

#s1 = random.uniform(0,1)
#s2 = random.uniform(0,1)
#s3 = random.uniform(0,1)


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

