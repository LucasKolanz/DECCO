##This file should go through the specified folders and extract the maximum and minimum aggregate radii.
##Radius of an aggregate is measured as the max distance from center of mass to particle


import os
import glob
import numpy as np
import subprocess
import check_for_errors as cfe
import check_for_errors_h5 as cfeh5
import h5py
import json
import sys

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u
import dataUtils as du

def findMonomerRadii(data_folder,relax=False):
	size = int(data_folder.split('/')[-3].split('_')[-1])
	try:
		pos,radii,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	except:
		print(f"no data in {data_folder}")
		return np.nan,np.nan
	if pos is None:
		return np.nan,np.nan


	return np.min(radii),np.max(radii)

def calcAggRadii(folder,relax=False):
	pos = du.getPos(folder,relax)
	if pos is None:
		return np.nan
	mass = du.getMass(folder,relax)
	COM = u.calcCOM(pos,mass)

	return max([np.linalg.norm(p-COM) for p in pos])

def main():


	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	curr_folder = os.getcwd() + '/'

	jobs = []

	# job = input_json["data_directory"] + job_folder + 'constant$a$/N_$n$/T_$t$/'
	# job_folder = 'jobsNovus/'
	# job = input_json["data_directory"] + job_folder + 'constrelax_$a$/N_$n$/T_$t$/'
	# jobs.append(job)
	job_folder = 'jobsCosine/'
	job = input_json["data_directory"] + job_folder + 'lognormrelax_$a$/N_$n$/T_$t$/'
	jobs.append(job)

	# print(job)
	# move_folder = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'


	attempts = [i for i in range(30)]
	# attempts = [1]

	N = [30,100,300]
	# N=[300]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	aggRadii = []
	minmonomerRadii = []
	maxmonomerRadii = []
	for job in jobs:
		relax = False
		if "relax" in job:
			relax = True
		print(f"relax: {relax}")
		for attempt in attempts:
			for n in N:
				for Temp in Temps:
					folder = job.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
					aggRadii.append(calcAggRadii(folder,relax))
					minradii,maxradii = findMonomerRadii(folder,relax)
					minmonomerRadii.append(minradii)
					maxmonomerRadii.append(maxradii)
	print(f"Max Agg radii: {np.nanmax(aggRadii)}")
	print(f"Min Agg radii: {np.nanmin(aggRadii)}")

	print(f"Max mono radii: {np.nanmax(maxmonomerRadii)}")
	print(f"Min mono radii: {np.nanmin(minmonomerRadii)}")

if __name__ == '__main__':
	main()