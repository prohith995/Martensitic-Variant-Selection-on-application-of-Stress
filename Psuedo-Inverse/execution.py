
import plots as p
from operator import itemgetter
import math
import other_functions as fn1

#strain = [[0.01,0,0],[0,0,0],[0,0,-0.01]] ## Rolling
strain = [[0.01,0,0],[0,-0.005,0],[0,0,-0.005]] ## Tension

n = 500

(poles, variants, energy) = fn1.solve(n,strain)

p.plot_pf(poles)
p.plot_vars(variants)
p.plot_energy(energy)
	