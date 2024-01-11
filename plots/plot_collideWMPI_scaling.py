import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from matplotlib.ticker import AutoMinorLocator

import sys
sys.path.append("/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab")
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
	base = os.getcwd() + "/jobs/initialScaling/"
	base = "/global/homes/l/lpkolanz/SpaceLab/jobs/"
	MPIbase = "/global/homes/l/lpkolanz/SpaceLab/testMPI/jobs/"
	MPIGPUbase = "/global/homes/l/lpkolanz/SpaceLab/testMPIGPU/jobs/"
	
	inds = np.arange(1,20)
	threads = [1,2,4,8,16,32,64,128]
	nodes = [1,2,4,8,16,32]
	thPerNode = 32
	folders = ["strongScaleCollide1"]#,"weakScaleGrowth1"]
	folders = ["strongScaleCollide_O2_1"]#,"weakScaleGrowth1"]
	folders = ["strongScaleCollide_O2_1200_1"]#,"weakScaleGrowth1"]
	folders = [base+"strongScaleCollide_O2_2400_1",MPIbase+"strongScaleCollideMPI2",MPIGPUbase+"strongScaleCollideMPIGPU1"]
	# inds = np.arange(1,3)

	times = np.zeros((len(folders),len(threads)),dtype=np.float64)
	times[:,:] = np.nan
	
	temp = 100
	for f_i,folder in enumerate(folders):
		if folder == folders[0]:
			th = threads
			timeBase = folder + "/thread_"
		else:
			th = nodes
			timeBase = folder + "/node_"
		# lowest_index = 100
		th = threads
		#if folder == folders[1]:
		#	th = th[:-1]
		for t_i,t in enumerate(th):
			timeFile = timeBase + "{}/time.csv".format(t,t)
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
				print(str(time) + " " + str(f_i) + " " + str(t_i))
			except FileNotFoundError:
				continue

	# times[0] = 
	# times[1] = 
	# times = np.copy(times)
	# times[0,:] = times[0,0]/times[0,:]
	# times[1,:] = times[1,0]/times[1,:]

	print(times)	
	# print(threads)		

	total_x_data = []
	total_x_data.extend(threads)
	total_x_data.extend([node*thPerNode for node in nodes])
	total_x_data = np.unique(total_x_data)
	# print(nodes)
	# print([node*thPerNode for node in nodes])
	print(total_x_data)

	b = times[0,0]*threads[0]

	# ideal = np.power(2,-1*np.log(threads)+b)
	ideal = b/total_x_data
	# print(ideal)

	# efficiency = (times[0,6])/(times[0,2]/threads[2])
	# print(efficiency)
	# print((times[0,5])/(times[0,2]/threads[2]))
	# print((times[0,3]/threads[3])/(times[0,2]/threads[2]))


	fig, ax = plt.subplots(1,1,figsize=(15,7))
	ax.tick_params(axis='x',labelsize=fontsize)
	ax.tick_params(axis='y',labelsize=fontsize)
	plt.rcParams.update({'font.size': fontsize})
	ax.set_xscale('log', base=2)
	ax.set_yscale('log', base=2)
	ax.minorticks_on()
	for f_i,folder in enumerate(folders):
		xdata = []
		title = "Strong Scaling of threads and MPI"
		if folder == folders[0]:
			typ = "Threads"
			c = 'green'
			m = '*'
			xdata = threads
		else:
			typ = "MPI"
			c = 'black'
			m = '.'
			xdata = [node*thPerNode for node in nodes]
			print("HERE")
			print(xdata)

		# print(inds)
		# print(times[f_i,:len(inds)])

		ax.plot(xdata,times[f_i,:len(xdata)],label='{} strong scaling'.format(typ),color=c,marker=m)
		
		# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
		# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
		# ax.grid()
		# ax.grid(visible=True,which = "minor")
	ax.plot(total_x_data,ideal,label="Ideal strong scaling")
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

	# ax.loglog(inds,times[f_i,:len(inds)])
	# ax.plot(inds,,label='multiCoreTest7')
	# ax.set_title(title)
	ax.set_xlabel("Number of Threads",fontsize=fontsize)
	ax.set_ylabel("Time (s)",fontsize=fontsize)
	ax.legend()
	plt.tight_layout()
	# plt.savefig("figures/{}loglogCollideScaling.png".format(title.split(' ')[0]))
	# plt.savefig("figures/{}loglogCollideScaling_O2.png".format(title.split(' ')[0]))
	plt.savefig("figures/MPIandThreadsStrongScaling.png")


	plt.show()

if __name__ == '__main__':
	main()
