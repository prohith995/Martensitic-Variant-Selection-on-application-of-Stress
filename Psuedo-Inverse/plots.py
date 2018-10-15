
import matplotlib.pyplot as plt

# Plotting pole figure
def plot_pf(poles):

	# Plooting the skeleton of the pole figure
	plt.rcParams.update({'font.size': 15})
	circle1 = plt.Circle((0, 0), 1, fill=False, linewidth = 1.8)
	fig = plt.gcf()
	ax = fig.gca()
	ax.set_aspect('equal')
	ax.add_artist(circle1)
	plt.title("0 0 1 Pole Figure for Tensile Strain - rolled texture")
	plt.plot((-1, 1), (0, 0), 'k-', linewidth = 1.8)
	plt.plot((0, 0), (-1, 1), 'k-', linewidth = 1.8)	

	# Plotting the poles
	x = [i[0] for i in poles]
	y = [i[1] for i in poles]

	plt.axis((-1,1,-1,1)) ## Axes limits
	plt.plot(x,y,".")
	plt.show()

# Plotting the bar chart of the number of times each variant occurs
def plot_vars(variants):
	n = [i for i in range(1,25)]
	count = [0]*24
	for i in n :
		count[i-1] = variants.count(i)
	plt.figure()
	plt.rcParams.update({'font.size': 25})
	plt.bar(n,count)
	plt.title("Preference of variants for tensile strain")
	#plt.ylim([0, 1800])
	plt.xlabel('Variant Number')
	plt.ylabel('Number of times occuring')
	plt.show()

# Plotting the histogram of the energy
def plot_energy(energy):
	plt.figure()
	plt.rcParams.update({'font.size': 25})
	plt.title("Histogram of strain energy for tension")
	plt.xlabel('Energy')
	plt.ylabel('Number of times occuring')
	#plt.ylim([0, 550])
	numBins = 100
	plt.hist(energy,numBins,color='green',alpha=0.8)
	plt.show()
