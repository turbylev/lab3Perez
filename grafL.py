import numpy as np
import matplotlib.pyplot as plt
try:	
	data=np.loadtxt ("Lag.txt")
	plt.plot(data[:,0], data[:,1])
	data=np.loadtxt ("Lag2.txt")
	plt.plot(data[:,0], data[:,1],'r.')
	plt.minorticks_on()
	plt.grid(color="black", which="major", linewidth=0.5)
	plt.grid(color="black", which="minor", linestyle=":", linewidth=0.5)
	plt.show()
except BaseException:
	print("error")
