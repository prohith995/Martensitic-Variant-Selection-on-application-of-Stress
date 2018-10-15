
import functions as fn
import other_functions as fn1
import plots as p
from operator import itemgetter
import math
import random
import numpy as np


strain = [[0.01,0,0],[0,0,0],[0,0,-0.01]]
#strain = [[0.01,0,0],[0,-0.005,0],[0,0,-0.005]] ## Tension

n = 5000

(poles, variants, energy) = fn1.simplex_run(n, strain)

p.plot_pf(poles)
p.plot_vars(variants)
p.plot_energy(energy)


## The following commented lines were used to verify the double dot product of stress and strain.
## Required changes are needed to be made in both the simplex.py and other functions.py

#table = fn1.simplex_obtain_stress(n, strain)
#strains = fn1.M(6)
#product = []
#for i in range(len(strains)) :
#	s1 = strains[i]
#	s2 = table[i+1]
#	s3 = fn.sixD_inverse(s1)
#	s4 = fn.sixD_inverse(s2)
#	product.append( [fn.double_dotproduct(s3,s4), i+1] )
#for i in range(len(product)) :
#	product[i] = [abs(product[i][0]), product[i][1]]
#product = sorted(product,key=itemgetter(0),reverse = True)
#for i in product :
#	print i

