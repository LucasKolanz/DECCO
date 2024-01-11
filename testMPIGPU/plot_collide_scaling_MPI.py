import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from matplotlib.ticker import AutoMinorLocator

import sys
sys.path.append("/global/homes/l/lpkolanz/SpaceLab/")
import utils as u
import porosity_FD as p
from matplotlib import ticker


# Function to generate minor ticks for log base 2 scale
def log_2_minor_ticks(min_val, max_val, num_per_major=4):
    min_pow = np.floor(np.log2(min_val))
    max_pow = np.ceil(np.log2(max_val))
    minor_ticks = []
    for i in range(int(min_pow), int(max_pow) + 1):
        major_tick = 2 ** i
        minor_spacing = major_tick / (num_per_major + 1)
        minor_ticks.extend([2 ** i + j * minor_spacing for j in range(1, num_per_major + 1)])
    return minor_ticks


def main():
	fontsize = 25
	# base = os.getcwd() + "/jobs/initialScaling/"
	base = "/global/homes/l/lpkolanz/SpaceLab/testMPIGPU/jobs/"
	base = "/pscratch/sd/l/lpkolanz/SpaceLab/testMPIGPU/jobs/"
	
	inds = np.arange(1,20)
	threads = [1,2,4,8,16,32,64,128]
	threads = [1,2,4,8,16,32]

	folders = ["strongScaleCollide1"]#,"weakScaleGrowth1"]
	folders = ["strongScaleCollide_O2_1"]#,"weakScaleGrowth1"]
	folders = ["strongScaleCollideSimple1"]#,"weakScaleGrowth1"]
	# folders = ["strongScaleCollideMPIonethread1"]#,"weakScaleGrowth1"]
	# inds = np.arange(1,3)

	times = np.zeros((len(folders),len(threads)),dtype=np.float64)
	times[:,:] = np.nan
	
	temp = 100
	for f_i,folder in enumerate(folders):
		# lowest_index = 100
		th = threads
		#if folder == folders[1]:
		#	th = th[:-1]
		for t_i,t in enumerate(th):
			timeFile = base + folder + "/node_{}/time.csv".format(t,t)
			print("========================================")
			print(timeFile)
			try:
				with open(timeFile,'r') as tF:
					lines = tF.readlines()
					print(lines)
				# if t == th[0]:
				# 	i = -1
				# 	while "ball,update time" in lines[i]:
				# 		i -= 1
				# 	lowest_index = lines[i].split(',')
				# 	print(lowest_index)
				# else:
				# 	continue
				try:
					time = float(lines[-1].split(',')[1][:-1])
				except:
					continue
				# lowest_index = int(lines[-1].split(',')[0])
				times[f_i,t_i] = time
			except FileNotFoundError:
				continue

	# times[0] = 
	# times[1] = 
	speedups = np.copy(times)
	# speedups[0,:] = speedups[0,0]/speedups[0,:]
	# speedups[1,:] = speedups[1,0]/speedups[1,:]

	print(speedups)	
	print(threads)		

	firstnonnan = 0
	b = speedups[0,firstnonnan]*threads[firstnonnan]

	# ideal = np.power(2,-1*np.log(threads)+b)
	ideal = b/threads
	print("ideal")
	print(ideal)
	# ideal[0:2] = np.nan

	# efficiency = (speedups[0,6])/(speedups[0,2]/threads[2])
	# print(efficiency)
	# print((speedups[0,5])/(speedups[0,2]/threads[2]))
	# print((speedups[0,3]/threads[3])/(speedups[0,2]/threads[2]))


	for f_i,folder in enumerate(folders):
		if folder == folders[0]:
			inds = threads
			title = "MPI Strong Scaling of sim_one_step"
		else:
			inds = threads[:-1]
			title = "MPI Weak Scaling of sim_one_step"
		# print(inds)
		# print(speedups[f_i,:len(inds)])
		fig, ax = plt.subplots(1,1,figsize=(15,7))
		ax.tick_params(axis='x',labelsize=fontsize)
		ax.tick_params(axis='y',labelsize=fontsize)
		plt.rcParams.update({'font.size': fontsize})
		ax.set_xscale('log', base=2)
		ax.set_yscale('log', base=2)
		ax.minorticks_on()

		ax.plot(inds,speedups[f_i,:len(inds)],label='Strong scaling',color='green',marker='*')
		ax.plot(inds,ideal,label="Ideal strong scaling")
		
		# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
		# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
		# ax.grid()
		# ax.grid(visible=True,which = "minor")
		ax.grid(visible=True,which = "major", linewidth = 1)
		ax.grid(visible=True,which = "minor", linewidth = 0.2)
		
		x_min, x_max = ax.get_xlim()
		y_min, y_max = ax.get_ylim()

		x_minor_ticks = log_2_minor_ticks(x_min, x_max,4)
		y_minor_ticks = log_2_minor_ticks(y_min, y_max,4)

		ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))
		ax.yaxis.set_minor_locator(ticker.FixedLocator(y_minor_ticks))

		ax.xaxis.set_minor_formatter(ticker.NullFormatter())
		ax.yaxis.set_minor_formatter(ticker.NullFormatter())
		# Set the custom formatter for the y-axis major ticks
		# ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_major_formatter))


		ax.grid(visible=True,which = "major", linewidth = 0.75,color='black')
		ax.grid(visible=True,which = "minor", linewidth = 0.2,color='black')

		# ax.loglog(inds,speedups[f_i,:len(inds)])
		# ax.plot(inds,,label='multiCoreTest7')
		ax.set_title(title)
		ax.set_xlabel("Number of Nodes (each w/ 1 GPU)",fontsize=fontsize)
		ax.set_ylabel("Time (s)",fontsize=fontsize)
		ax.legend()
		plt.tight_layout()
		# plt.savefig("figures/{}loglogCollideScaling.png".format(title.split(' ')[0]))
		# plt.savefig("figures/{}loglogCollideScaling_O2.png".format(title.split(' ')[0]))
		# plt.savefig("figures/{}MPIGPUloglogCollideScaling.png".format(title.split(' ')[0]))


	plt.show()

if __name__ == '__main__':
	main()
