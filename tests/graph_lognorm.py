import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv

def lognorm(a,sigma,a_max):
	return 1/(a*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(a/a_max)-sigma**2)**2/(2*sigma**2))

def main():
	data_lognorm = np.genfromtxt('lognorm_output.csv',delimiter=",",dtype=np.float64)
	
	xdata = np.linspace(1,5,100)
	ydata = lognorm(xdata,0.2,2)

	fig,ax = plt.subplots()

	ax.hist(data_lognorm,bins=100,density=True)
	ax.plot(xdata,ydata)
	plt.show()

if __name__ == '__main__':
	main()