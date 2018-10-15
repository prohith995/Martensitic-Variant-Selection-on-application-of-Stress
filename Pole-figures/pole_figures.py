
import numpy as np
import math
import matplotlib.pyplot as plt

from test1data import *
import macro_strain


# To transform the normal vectors
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

# Transforms from crystal frame to sample frame
def transform_inverse(Q,tensor):
	t1 = np.dot(Q.transpose(),tensor)
	tf = np.dot(t1, Q)
	return tf

# Returns the poles (in sample frame) taking the euler angles of the polycrystal as input
def run(pref_var,euler_angles):

	# Obtaining the normal directions of all the variants in the preferred variants list
	norm_pref_var1 = [p1[i[0]-1] for i in pref_var] 	
	norm_pref_var = []
	for i in norm_pref_var1:
		norm_pref_var.append(i)
		norm_pref_var.append([-i[0],-i[1],-i[2]])
	
	Q = np.asarray(transformation(euler_angles))
	# Obtaining the normal vectors in the sample frame
	norm_pref_var_samp_fram = [transform_inverse(Q,i) for i in norm_pref_var]
	# print norm_pref_var_samp_fram

	poles = []
	for i in norm_pref_var_samp_fram :
		if i[2] >= 0 : ## To see if the pole exists in the pole figure. i.e. it is in the northen hemisphere
			x = i[0]/(i[2]+1)
			y = i[1]/(i[2]+1)
			poles.append([x,y])
	#print poles
	return poles
	
# Euler angles of polycrystal
f = open("rand5000.txt",'r')
n = 50
poles = []
for i in range(n):
	a = f.readline().split()
	euler_angles = [float(a[0]), float(a[1]), float(a[2])]
	#euler_angles = [128.99, 48.27, 24.27]
	#euler_angles = [0, 45, 0]
	print "\n\n", i
	pref_var = macro_strain.obtain_preferred_variants(euler_angles)
	for i in pref_var :
		print i[0], 
	#print pref_var
	poles.extend(run(pref_var,euler_angles))
	#print poles

## Plotting pole figure

## Plotting the skeleton of the pole figure
circle1 = plt.Circle((0, 0), 1, fill=False)
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect('equal')
ax.add_artist(circle1)
plt.plot((-1, 1), (0, 0), 'k-')
plt.plot((0, 0), (-1, 1), 'k-')	

# Plotting the poles
#print poles
x = [i[0] for i in poles]
y = [i[1] for i in poles]

#print x
#print y

plt.axis((-2,2,-2,2)) ## Axes limits
plt.plot(x,y,"o")

plt.show()




