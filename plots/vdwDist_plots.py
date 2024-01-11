import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []

	attempt = [0,1,2,3,4,5]
	attempt = [6,7,8]
	attempt = [0,1,2,3,4,5,6,7,8]
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

	data = []
	# for i in range(0,len(path),3):
	for i in range(0,len(path)):
		# fig, ax = plt.subplots(1,len(attempt),figsize=(10,7))
		# fig, ax = plt.subplots(1,1,figsize=(10,7))
		# ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
		# if len(attempt) > 1:
		# 	ax = ax.flatten()
		# else:
		# ax = [ax]
		# for i in range(j,j+len(attempt)):
			
		# distfile = path[i] + fileprefix[i] + "distData.csv"
		# vdwfile = path[i] + fileprefix[i] + "vdwData.csv"
		inputfile = path[i] + "input.json"
		
		# distdata = np.loadtxt(distfile,delimiter=',',dtype=np.float64,skiprows=1)
		# vdwdata = np.loadtxt(vdwfile,delimiter=',',dtype=np.float64,skiprows=1)[:,:3]
		# dist = distdata
		# vdw = [np.sqrt(v[0]**2 + v[1]**2 + v[2]**2) for v in vdwdata[:-1]]
		
		inputs = json.load(open(inputfile))
		h_min = float(inputs['h_min'])
		r1 = float(inputs['note'].split(',')[0].split('=')[1])
		r2 = float(inputs['note'].split(',')[1].split('=')[1])
		scaleBalls = float(inputs['scaleBalls'])
		Ha = float(inputs['Ha'])

		dist = np.linspace(-1e-5,4e-5,num=int(1e5))
		distCalc = np.where(dist > h_min*scaleBalls, dist, h_min*scaleBalls)
		vdwForce = (Ha*(64/6))*(r1**3*r2**3*(distCalc+r1+r2)/ \
					((distCalc**2 + 2*r1*distCalc + 2*r2*distCalc)**2*\
					(distCalc**2 + 2*r1*distCalc + 2*r2*distCalc + 4*r1*r2)**2))
		data.append(vdwForce)
		# print(dist[0:5])

		# fig = plt.figure()
		# ax = plt.axes()

		# # print(end)
		# # print(dist[start:end])
		# # print(dist)
		# # print(vdw[start:end])

		# ax.plot(dist,vdwForce,color='g',label="r1={}, r2={}, hmin={}".format(r1,r2,h_min*scaleBalls))
		# ax.legend()
		# ax.set_xlabel("dist [cm]")
		# ax.set_ylabel("VDW force [ergs]")
		# ax.set_title("Theoretical VDW vs distance between balls")
		# # ax[i].set_xlim((np.min(dist),2e-5))
		# # ax[i].set_xlim((1.9674e-5,1.96740000000000000001e-5))
		# ax.axvline(h_min*scaleBalls)

		# plt.tight_layout()

	for i in range(3):
		fig = plt.figure()
		ax = plt.axes()
		fig1 = plt.figure()
		ax1 = plt.axes()
		
		ax.plot(dist,data[len(data)-3+i]/data[i],color='g',label="data[{}]/data[{}]".format(len(data)-3+i,i))
		ax.legend()
		ax.set_xlabel("dist [cm]")
		ax.set_ylabel("VDW force ratio")
		ax.set_title("Theoretical VDW ratio vs distance between balls")
		# ax[i].set_xlim((np.min(dist),2e-5))
		# ax[i].set_xlim((1.9674e-5,1.96740000000000000001e-5))
		ax.axvline(h_min*scaleBalls)
		ax.set_yscale("log")

		ax1.plot(dist,data[len(data)-3+i]/data[i+3],color='g',label="data[{}]/data[{}]".format(len(data)-3+i,i+3))
		ax1.legend()
		ax1.set_xlabel("dist [cm]")
		ax1.set_ylabel("VDW force ratio")
		ax1.set_title("Theoretical VDW ratio vs distance between balls")
		# ax[i].set_xlim((np.min(dist),2e-5))
		# ax[i].set_xlim((1.9674e-5,1.96740000000000000001e-5))
		ax1.axvline(h_min*scaleBalls)
		ax1.set_yscale("log")

		plt.tight_layout()
		
	plt.show()

if __name__ == '__main__':
	main()