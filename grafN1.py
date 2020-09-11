import numpy as np
import matplotlib.pyplot as plt
try:	
	data=np.loadtxt ("New1.txt")
	plt.plot(data[:,0], data[:,1], marker='.')
	plt.gcf().canvas.set_window_title("Polinom1 N'utona")
	plt.minorticks_on()
	plt.grid(color="k", which="major", linewidth=1)
	plt.grid(color="k", which="minor", linestyle=":", linewidth=0.5)
	plt.show()
except BaseException:
	print("error")
