import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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


# Set the custom formatter for the y-axis major ticks
def y_major_formatter(x, pos):
    power = int(np.log2(x))
    return f'$2^{{{power}}}$' if power % 2 != 0 else ''


def main():
	fontsize = 25
	base = os.getcwd() + "/jobs/initialScaling/"
	base = "/global/homes/l/lpkolanz/SpaceLab/testMPI/jobs/"
	

	inds = np.arange(1,20)
	nodes = [1,4,16,64]
	folders = []#["weakScaleCollideMPI1/thread_4/","weakScaleCollideMPI1/thread_16/","weakScaleCollideMPI1/thread_64/"]
	for node in nodes:
		folders.append("weakScaleCollideMPI1/node_{}/".format(node))
	# inds = np.arange(1,3)

	times = np.zeros((1,len(nodes)),dtype=np.float64)
	times[:,:] = np.nan

	
	temp = 100
	for f_i,folder in enumerate(folders):

		timeFile = base + folder + "time.csv"
		print("========================================")
		print(timeFile)
		try:
			with open(timeFile,'r') as tF:
				lines = tF.readlines()
			# if t == th[0]:
			# 	i = -1
			# 	while "ball,update time" in lines[i]:
			# 		i -= 1
			# 	lowest_index = lines[i].split(',')
			# 	print(lowest_index)
			# else:
			# 	continue
			print(lines)
			try:
				time = float(lines[-1].split(',')[1][:-1])
			except:
				continue
			# lowest_index = int(lines[-1].split(',')[0])
			times[0,f_i] = time
		except FileNotFoundError:
			continue

	# times[0] = 
	# times[1] = 
	speedups = np.copy(times)
	# speedups[0,:] = speedups[0,0]/speedups[0,:]
	# speedups[1,:] = speedups[1,0]/speedups[1,:]

	print(speedups)	
	print(speedups[0,0:2])	
	print(nodes)		

	b = speedups[0,0]

	# ideal = np.power(2,-1*np.log(nodes)+b)
	ideal = [b for i in nodes]
	# ideal = b/nodes
	print(ideal)
	# ideal[0:2] = np.nan

	efficiency = (speedups[0,2])/(speedups[0,0]/nodes[0])
	print(efficiency)


	# for f_i,folder in enumerate(folders):
		# inds = nodes
	title = "MPI Weak Scaling of sim_one_step"

	fig, ax = plt.subplots(1,1,figsize=(15,7))
	plt.rcParams.update({'font.size': fontsize})
	ax.tick_params(axis='x',labelsize=fontsize)
	ax.tick_params(axis='y',labelsize=fontsize)
	ax.set_xscale('log', base=2.0)
	ax.set_yscale('log', base=2.0)
	ax.minorticks_on()
	
	ax.bar(nodes[0],speedups[0,0],label='Weak scaling',color='g',width=1)
	ax.bar(nodes[1],speedups[0,1],color='g',width=4)
	ax.bar(nodes[2],speedups[0,2],color='g',width=16)
	# ax.bar(nodes[2],speedups[0,2],color='g',width=64)
	# ax.plot(nodes,speedups[0,:],label='Weak scaling',marker='.')
	ax.plot(nodes,ideal,label="Ideal weak scaling",linewidth=2)
	ax.set_ylim((1,2**6))
	# ax.set_xlim((2**1.5,2**6.6))

	ax.tick_params(axis='x',labelsize=fontsize)
	ax.tick_params(axis='y',labelsize=fontsize)
	plt.rcParams.update({'font.size': fontsize})

	y_major_ticks = [2**i for i in range(int(np.log2(1)), int(np.log2(np.max(speedups[0,0:2]))) + 1)]
	ax.set_yticks(y_major_ticks)

	x_min, x_max = ax.get_xlim()
	y_min, y_max = ax.get_ylim()

	x_minor_ticks = log_2_minor_ticks(x_min, x_max,4)
	y_minor_ticks = log_2_minor_ticks(y_min, y_max,4)

	ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))
	ax.yaxis.set_minor_locator(ticker.FixedLocator(y_minor_ticks))

	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	ax.yaxis.set_minor_formatter(ticker.NullFormatter())
	# Set the custom formatter for the y-axis major ticks
	ax.yaxis.set_major_formatter(plt.FuncFormatter(y_major_formatter))
	# ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_major_formatter))

	ax.grid(visible=True,which = "major", linewidth = 0.75,color='black')
	ax.grid(visible=True,which = "minor", linewidth = 0.2,color='black')
	# ax.loglog(inds,speedups[f_i,:len(inds)])
	# ax.plot(inds,,label='multiCoreTest7')
	ax.set_title(title)
	ax.set_xlabel("Number of Nodes",fontsize=fontsize)
	ax.set_ylabel("Time (s)",fontsize=fontsize)
	ax.legend(loc='lower right')
	plt.tight_layout()
	plt.savefig("figures/{}MPIloglogCollideScaling.png".format(title.split(' ')[0]))


	plt.show()

if __name__ == '__main__':
	main()
