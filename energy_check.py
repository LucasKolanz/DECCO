"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file verifies energy conservation in both the energy output file and simData file for given folders.

"""





import sys
import glob
import os
import json
import numpy as np
import matplotlib.pyplot as plt

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

def calc_KE(vel,w,mass,moi):
	#tranlational KE
	KE = 0.5 * np.einsum('i,i->', mass, np.einsum('ij,ij->i', vel, vel))
	#rotational KE
	KE += 0.5 * np.einsum('i,i->', moi, np.einsum('ij,ij->i', w, w))
	return KE

def calc_transKE(vel,mass):
	#tranlational KE
	KE = 0.5 * np.einsum('i,i->', mass, np.einsum('ij,ij->i', vel, vel))
	return KE

def calc_rotKE(w,moi):
	#rotational KE
	KE = 0.5 * np.einsum('i,i->', moi, np.einsum('ij,ij->i', w, w))
	return KE

#calculates VDW and spring PE given the position of each ball pos, the radius of each ball radius, 
#the Hamakers constant between pairs HA, and spring constant K
#NOT FINISHED
def calc_PE(pos,radius,HA,K):
	balls = len(radius)
	PE = 0
	#loop over all particle pairs
	for i in range(1,balls):
		for j in range(0,i):
			dist = np.linalg.norm(pos[i] - pos[j])
			overlap = radius[i] + radius[j] - dist

			#overlapping so use spring PE
			if overlap > 0:
				PE += 0.5*K*overlap**2
			#not overlapping so use VDW PE
			else:
				pass

	exit(0)

def verify_simData(directory,HA=0):
	indices = u.get_all_indices(directory,checkpoint=True)
	indices = [4]
	relax = ("relax" in directory)
	for i in indices:
		time,_,_,_,_,_ = u.get_energy(directory,data_index=i,relax=relax)
		radius,mass,moi = u.get_constants(directory,i,relax=relax)
		pos,vel,w = u.get_simData(directory,i,relax)
		timesteps = pos.shape[0]

		start = 100
		stop = -1

		rotKE = np.zeros((timesteps),dtype=np.float64)
		transKE = np.zeros((timesteps),dtype=np.float64)

		# PE = np.zeros((timesteps),dtype=np.float64)

		for step in range(timesteps):
			transKE[step] = calc_transKE(vel[step],mass)
			rotKE[step] = calc_rotKE(w[step],moi)

			# PE[step] = calc_PE(pos[step],radius,HA=0.1,K=0.1)
		plt.figure(figsize=(10, 6))
		# plt.plot(time, KE, label='Kinetic Energy')
		plt.plot(time[start:stop], rotKE[start:stop], label='Rotational KE')
		plt.plot(time[start:stop], transKE[start:stop], label='Translational KE')
		plt.plot(time[start:stop], transKE[start:stop]+rotKE[start:stop], label='Total KE')
		# plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
		plt.ylim(0,2.2e-15)
		plt.xlabel('Time (s)')
		plt.ylabel('Energy (ergs)')
		# plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
		plt.legend()
		plt.grid(True)
		# plt.show()

		plt.figure(figsize=(10, 6))
		# plt.plot(time, KE, label='Kinetic Energy')
		plt.plot(time[start:stop], w[:,0,2][start:stop], label='w0')
		plt.plot(time[start:stop], w[:,1,2][start:stop], label='w1')
		# plt.plot(time, w[:,2,2], label='w2')
		# plt.plot(time, transKE, label='Translational KE')
		# plt.plot(time, transKE+rotKE, label='Total KE')
		# plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
		plt.xlabel('Time (s)')
		plt.ylabel('angular velocity (rad/s)')
		# plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
		plt.legend()
		plt.grid(True)


def verify_energy(directory):
	indices = u.get_all_indices(directory,checkpoint=True)
	indices = [4]
	relax = ("relax" in directory)
	for i in indices:
		time,PE,KE,E,p,L = u.get_energy(directory,data_index=i,relax=relax)
		
		start = 1
		stop = -1

		# Plot rotational kinetic energy versus time
		plt.figure(figsize=(10, 6))
		# plt.plot(time[start:stop], KE[start:stop], label='Kinetic Energy')
		plt.plot(time[start:stop], PE[start:stop], label='Potential Energy')
		plt.plot(time[start:stop], E[start:stop], label='Total Energy')
		# plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
		plt.xlabel('Time (s)')
		plt.ylabel('Energy (ergs)')
		# plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
		plt.legend()
		plt.grid(True)

		# print(KE[1])
		# print(PE[1])
		# print(E[1])

		# plt.figure(figsize=(10, 6))
		# # plt.plot(time, KE, label='Kinetic Energy')
		# # plt.plot(time, L, label='L')
		# plt.plot(time, p, label='p')
		# # plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
		# plt.xlabel('Time (s)')
		# plt.ylabel('Energy (ergs)')
		# # plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
		# plt.legend()
		# plt.grid(True)
		# # plt.show()



def verify_contactRadius(directory):
	data = np.loadtxt(directory+"contactRadii.txt",delimiter=",",dtype=np.float64)

	a = data[:,0]
	overlap = data[:,1]
	force = data[:,2]

	# print(a)

	start = 1
	stop = 1000

	# Plot rotational kinetic energy versus time
	plt.figure(figsize=(10, 6))
	plt.plot(list(range(len(a[start:stop]))), a[start:stop], label='contact radius')
	plt.plot(list(range(len(a[start:stop]))), overlap[start:stop], label='overlap')

	# plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
	plt.xlabel('step')
	plt.ylabel('a (cm)')
	# plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
	plt.legend()
	plt.grid(True)


def verify_displacements(directory):
	data = np.loadtxt(directory+"out.txt",delimiter=",",dtype=np.float64)

	slidingDisp = data[:,0]
	rollingDisp = data[:,1]
	overlap = data[:,2]

	# print(a)

	start = 0
	stop = 10000

	# Plot rotational kinetic energy versus time
	plt.figure(figsize=(10, 6))
	plt.plot(list(range(len(slidingDisp[start:stop]))), slidingDisp[start:stop], label='sliding')
	plt.plot(list(range(len(slidingDisp[start:stop]))), rollingDisp[start:stop], label='rolling')
	# plt.plot(list(range(len(slidingDisp[start:stop]))), overlap[start:stop], label='overlap')

	# plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
	plt.xlabel('step')
	plt.ylabel('a (cm)')
	# plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
	plt.legend()
	plt.grid(True)




if __name__ == '__main__':
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	directory = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_branch/SpaceLab_data/jobs/JKRBPCA0/N_300/T_3/"
	directory = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_branch/SpaceLab_data/jobs/JKRTest_5/"
	relax = ("relax" in directory)

	verify_energy(directory)
	verify_simData(directory)
	# verify_contactRadius(directory)
	# verify_displacements(directory)
	plt.show()

	
