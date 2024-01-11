import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# initialization function: plot the background of each frame
# def init(line):
# 	line.set_data([], [])
# 	return line,

def update_lines(num, data, line):
	line.set_data(data[0,:num],data[1,:num])
	return line

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []

	attempt = [3,4,5]
	attempt = [5]
	attempt = [8]
	# attempt = [0]
	ref_file = ""
	for i in attempt:
		# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/cuttoff_test/c_{}/".format(i))
		path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/vdwDist{}/".format(i))
		for file in os.listdir(path[-1]):
			if len(file) > len(ref_file) and file[-4:] == ".csv":
				ref_file = file

		ref_file = "_".join(ref_file.split('_')[:-1]) + "_"
		fileprefix.append(ref_file)

	# print(data)
	# print(data[1])
	# # exit(0)
	new_data = False
	animate = False


	# for i in range(0,len(path),3):
	for i in range(0,len(path)):
		sav_ani = path[i] + 'vdwAni.csv'

		has_ani = os.path.exists(sav_ani)

		inputfile = path[i] + "input.json"
		inputs = json.load(open(inputfile))
		h_min = float(inputs['h_min'])
		r1 = float(inputs['note'].split(',')[0].split('=')[1])
		r2 = float(inputs['note'].split(',')[1].split('=')[1])
		scaleBalls = float(inputs['scaleBalls'])
		# print("For i={}, h_min={}".format(i,h_min))
		
		if new_data or not has_ani:
			distfile = path[i] + fileprefix[i] + "distData.csv"
			vdwfile = path[i] + fileprefix[i] + "vdwData.csv"
			
			distdata = np.loadtxt(distfile,delimiter=',',dtype=np.float64,skiprows=1)
			vdwdata = np.loadtxt(vdwfile,delimiter=',',dtype=np.float64,skiprows=1)[1:,-3:]
			if len(distdata.shape) > 1:
				dist = distdata[:,-1]
			else:
				dist = distdata

			vdw = np.linalg.norm(vdwdata,axis=1)

			# start = np.where(dist>=0.0)[0][0]
			# end = np.where(dist>=2e-5)[0][0]

			start = 0
			end = -1

			xdata = dist[start:end]-(r1+r2)
			ydata = vdw[start:end]
			# print(len(xdata))
			# print(len(ydata))
			# exit(0)
			data = np.array([xdata,ydata])
			np.savetxt(sav_ani,data,delimiter=',')
		else:
			data = np.loadtxt(sav_ani,delimiter=',',dtype=np.float64)
			if len(data) > 1e6:
				data = data[::int(len(data)/1e3)]
			xdata = data[0,:]
			ydata = data[1,:]

		fig = plt.figure()
		minx = min(xdata)
		maxx = max(xdata)
		miny = min(ydata)
		maxy = max(ydata)
		diffx = maxx-minx
		diffy = maxy-miny
		# ax = plt.axes(xlim=(-0.4e-6,-.2e-6),ylim=(miny-diffy/100,maxy+diffy/100))
		# ax = plt.axes(xlim=(-3.261e-7,-3.2475e-7),ylim=(2,2.25))


		if animate:
			ax = plt.axes(xlim=(minx-diffx/100,maxx+diffx/100),ylim=(miny-diffy/100,maxy+diffy/100))
			ax.set_title('r1={}, r2={}, hmin={}'.format(r1,r2,h_min*scaleBalls))
			ax.set_xlabel("dist [cm]")
			ax.set_ylabel("VDW force [ergs]")
			ax.axvline(h_min*scaleBalls,label='h_min')
			frames = np.arange(0,len(xdata),int(len(xdata)/1000))
			frames = np.append(frames,len(xdata))

			frames = frames[-5:]
			line = ax.plot(data[0,0:1],data[1,0:1],label="VDW Force",color='g',\
						markevery=[-1],marker='*')

			line_ani = animation.FuncAnimation(fig=fig, func=update_lines, init_func=None,\
								frames=frames, fargs=(data, line[0]),blit=False,interval=10)
			# fig.canvas.draw()
			line_ani.repeat = False
		else:
			ax = plt.axes()
			ax.set_title('r1={}, r2={}, hmin={}'.format(r1,r2,h_min*scaleBalls))
			ax.set_xlabel("dist [cm]")
			ax.set_ylabel("VDW force [ergs]")
			ax.axvline(h_min*scaleBalls,label='h_min')
			line = ax.plot(data[0],data[1],label="VDW Force",color='g',\
						markevery=[-1],marker='*')
		ax.set_yscale("log")
		ax.legend()
		plt.tight_layout()

	plt.show()

if __name__ == '__main__':
	main()