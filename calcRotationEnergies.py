"""
This file was originally written for SpaceLab/DECCO to check the seeds of all jobs in folders matching a specified pattern.

Author: Lucas Kolanz

This file is meant to go through all jobs contained in a specified file pattern and print out their seeds. All of these seeds should
be different, and this file is meant to help verify that.

TODO: use the seeds from input.json, not seedFile.txt. Save all seeds in a list to verify they are all unique.

"""





import os
import json
import glob
import numpy as np
import subprocess
import sys

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u



def main():

	#Open SpaceLab default file for directory information
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_template = input_json["data_directory"] + 'jobs/' + 'BAPA' + '{a}/N_{n}/T_{t}/'

	attempts = [i for i in range(10)]

	# N = [30,100,300]
	N=[300]

	# Temps = [3,10,30,100,300,1000]
	Temps = [1000]

	k = 1.380649/10**16  #boltzmann constant in g*cm^2/(s^2*K)

	for a in attempts:
		for n in N:
			for t in Temps:
				folder = job_template.replace("{a}",str(a)).replace("{n}",str(n)).replace("{t}",str(t))
				files = u.find_files(folder,"*simData.csv")	
				ke = []
				masses = u.get_masses(folder)
				for file in files:
					print(file)
					pos = u.get_all_pos_data(file)
					vel = u.get_all_vel_data(file)
					# radii = u.get_radii(folder)
					# u.plot_3d_spheres(data[-1],radii)
					# exit(0)
					# pos = np.ones(shape = (5,10,3))
					# vel = np.ones(shape = (5,10,3))
					# # masses = np.ones(shape = (10))
					# pos = np.array([[[1,0,0],[-1,0,0],[0,0,0]],[[1,0,0],[-1,0,0],[0,0,0]]])
					# vel = np.array([[[0,1,0],[0,-1,0],[0,0,0]],[[0,1,0],[0,-1,0],[0,0,0]]])
					# masses = np.array([1,1,1])

					ke.extend(u.calc_rotational_kinetic_energy(pos, vel, masses))

				u.plot_rotational_kinetic_energy(np.array(ke), 1.5 * k * t)


if __name__ == '__main__':
	main()