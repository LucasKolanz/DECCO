import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab")
import utils as u
import porosity_FD as p


def fixLabel(label):
	if len(label.split(',')) == 2:
		return label 
	prevnum = -1
	firstnum = 0
	count = 0
	newlabel = label[0] 
	labelSplit = label.split(',')
	for i,numstr in enumerate(labelSplit):
		num = int(numstr)
		
		if num == prevnum + 1:
			if numstr != labelSplit[0]:
				# if numstr != labelSplit[-1]:
				if numstr == labelSplit[-1] or int(labelSplit[i+1]) != num+1:
					newlabel += '-' + str(num)
				
			count += 1
		else:
			newlabel += ',' + str(num)

		prevnum = num

	return newlabel

def main():
	base = os.getcwd() + "/jobs/"

	inds = np.arange(1,20)
	# threads = [1,2,4,8,16,32,64]
	# folders = ["strongScaleGrowth1/thread_1/","affinityTests_th2_1/","affinityTests_th4_1/","affinityTests_th8_1/"]
	folders = ["strongScaleGrowth_O2_1/thread_1/","affinityTests_th2_O2_1/","affinityTests_th4_O2_1/","affinityTests_th8_O2_1/"]
	# inds = np.arange(1,3)

	times = [] 
	xlabels = [] 

	temp = 100
	for f_i,folder in enumerate(folders):
		# lowest_index = 100
		if f_i == 0:
			xlabels.append("1 th")
			try:
				with open(base+folder+"/time.csv",'r') as tF:
					lines = tF.readlines()
			except FileNotFoundError:
				times.append(np.nan)
				continue
			try:
				time = float(lines[1].split(',')[1][:-1])
			except:
				time = np.nan
			times.append(time)
		else:
			# print("HERE: " + folder)
			for subdir, dirs, files in os.walk(base+folder):
				for directory in dirs:
					threads = folder.split('_')[1][2:]
					affinity = directory.split('_')[1].replace('-',',')
					xlabels.append(fixLabel(affinity))
					try:
						with open(base+folder+directory+"/time.csv",'r') as tF:
							lines = tF.readlines()
					except FileNotFoundError:
						times.append(np.nan)
						continue
					try:
						time = float(lines[-1].split(',')[1][:-1])
					except:
						time = np.nan
					times.append(time)

	print(xlabels)


	title = "Affinity of sim_one_step"
	# print(inds)
	# print(speedups[f_i,:len(inds)])
	fig, ax = plt.subplots(1,1,figsize=(15,7))
	ax.plot(range(0,len(xlabels)),times,marker='*')
	ax.set_xticks(range(0,len(xlabels)))
	ax.set_xticklabels(xlabels)
	# ax.plot(inds,,label='multiCoreTest7')
	ax.set_title(title)
	ax.set_xlabel("Affinity")
	ax.set_ylabel("sim_one_step time (s)")
	# ax.legend()
	# plt.savefig("figures/affinity.png")
	plt.savefig("figures/affinity_O2.png")


	# plt.show()

if __name__ == '__main__':
	main()
