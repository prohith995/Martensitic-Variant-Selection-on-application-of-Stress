## Date : 23/02/2017
## This program was written by Rohith.P as part of BTP
## The class takes the objective function and constraints and gives out the outputs as required (Outputs have to be customized)
## Special cases or Exceptions have not yet been added

import random

class table :

	# Initiating the table
	def __init__(self, obj):
		self.rows = []
		self.rows.append([obj,0]) # Objective Function

	# Constraint of the type <=
	def constraint (self,lhs, rhs):
		self.rows.append([lhs,rhs])

	# This function adds the slack identity matrix to the simplex table
	def add_slack (self):
		n = len(self.rows)-1
		self.rows[0][0].extend([0]*n)
		for i in range(n):
			add = [0]*n
			add[i] = 1
			self.rows[i+1][0].extend(add)

	# This function performs one row operations required for one iteration of the simplex, with the pivot indices as input
	def row_op (self, pivot):
		
		# Pivot details
		p = pivot[0]
		q = pivot[1]
		pvt = pivot[2]

		# Performs the pivot row's row operation
		for j in range(len(self.rows[p][0])):
			self.rows[p][0][j] = self.rows[p][0][j]*1.0 / pvt
		self.rows[p][1] = self.rows[p][1]*1.0 / pvt
		
		# Performs row operation for remaining rows
		for i in range (len(self.rows)):
			if (i!=p) :
				c = self.rows[i][0][q]	# The constant required fro each row opn
				for j in range(len(self.rows[i][0])):
					self.rows[i][0][j] = self.rows[i][0][j] - c*self.rows[p][0][j]
				self.rows[i][1] = self.rows[i][1] - c*self.rows[p][1]

	# This function calculates the pivot
	def pivot (self):

		# Pivot column
		a = min(self.rows[0][0])
		pc = self.rows[0][0].index(a)
		pr = None
		if (a<0) :			
			# Pivot row calculation, minimum ratio test
			mini = -1
			for i in range(len(self.rows)-1,0,-1) :
				if (self.rows[i][0][pc] >= 0.0000001):
					ratio = self.rows[i][1]/self.rows[i][0][pc]
					if mini < 0 :
						if ratio >= 0 :
							pr = i
							mini = ratio
					else :
						if ratio <= mini and ratio >= 0 :
							pr = i
							mini = ratio

			if pr != None :
				return ([pr,pc,self.rows[pr][0][pc]])
		
			return None

	# Runs the simplex algortihm till optimum is found
	def run (self,constraints):
		for i in constraints :
			self.constraint(i[0],i[1])
		self.add_slack()

		a = self.pivot()
		iteration = 0
		while (a!=None): # If None, Optimal
			zz = self.rows[0][1]
			self.row_op(a)

			a = self.pivot()
			iteration = iteration + 1

		variants = self.out_positive()
		print self.rows[0][1]
		return variants, self.rows[0][1]
#		return self.out_cols(0,12) ## For the stress verification

	# Prints the present simplex table
	def out(self):
		for i in self.rows :
			print i

	# Prints and returns the active variants
		for i in range(len(self.rows[0][0])-len(self.rows)+1 , len(self.rows[0][0])) :
			if self.rows[0][0][i] > 0 :
				energy = energy + self.rows[0][0][i]
				print self.rows[0][0][i] , i-5
				variants.append(i-5)
		return variants