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
import check_for_errors as cfe

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def transfer(source,dest):
	if not os.path.exists(dest):
		os.makedirs(dest)
	os.system(f"cp {source}* {dest}.")

def get_attempt(folder):

	folder_base = folder 
	i = 1000
	folder = folder_base.replace("$a$",str(i))
	# print(folder+"timing.txt")
	# print(os.path.exists(folder+"timing.txt"))
	while os.path.exists(folder+"timing.txt"):
		i += 1
		folder = folder_base.replace("$a$",str(i))
		# print(folder+"timing.txt")
		# print(os.path.exists(folder+"timing.txt"))

	# print(folder)
	return i


def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
	data_dir = input_json["data_directory"]

	job_folder = data_dir + 'jobsNovus/'##FOR LOCAL
	job_folder = data_dir + 'jobsCosine/'##FOR LOCAL

	base_job = job_folder + 'const$a$/N_$n$/T_$t$/'
	base_job = job_folder + 'lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[30]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	for n in N:
		for t in Temps:
			print(f"====================  n={n} and t={t}  ====================")
			for a in attempts:
				folder = base_job.replace("$a$",str(a)).replace("$n$",str(n)).replace("$t$",str(t))

				if os.path.exists(folder+"seedFile.txt"):
					with open(folder+"seedFile.txt",'r') as f:
						print(f.read())


if __name__ == '__main__':
	main()