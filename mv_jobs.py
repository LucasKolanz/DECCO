##restart jobs that have errored

#Error 1: restart failed, number of balls didn't carry over so there are an incorrect 
#			number of colums in at least 1 sim_data output
#Error 2: sim fails to write all the data (so there aren't the correct number of colums in simData)
#Error *: possible signed integer over flow in number of steps for a sim
#			Indicated by a specific sequence at the end of sim_errors.txt
#Error general: did we get "Simulation complete!" within the last 10 lines of sim_error.log

import os
import glob
import numpy as np
import subprocess
import check_for_errors as cfe

# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def check_error(job_base,error,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	errors = []
	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts
	# N=[30]
	# Temps = [3]
	valid_count = 0

	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				if os.path.exists(job):
					if os.path.exists(job+"timing.txt"):
						valid_count += 1
					output = error(job)
					print(f"error: {output}")
					if not output:
						errors.append(job)
				else:
					print(f"Job doesn't exist: {job}")

	# print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return errors


def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	curr_folder = os.getcwd() + '/'


	sink_job_folder = 'jobs/'###FOR COSINE
	source_job_folder = 'erroredJobs/'

	sink_job = input_json["data_directory"] + sink_job_folder + 'lognorm$a$/N_$n$/T_$t$/'
	source_job = input_json["data_directory"] + source_job_folder + 'lognorm$a$/N_$n$/T_$t$/'
	# move_folder = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	attempts = [18]

	N = [30,100,300]
	N=[30]

	Temps = [3,10,30,100,300,1000]
	Temps = [3]

	error_folders = check_error(source_job,cfe.error4,N,Temps,attempts)

	if len(error_folders) == 0:
		print("No errors found")

	for folder in error_folders:
		command = f"mv {folder}* {folder.replace(source_job_folder,sink_job_folder)}."
		print(command)
		# os.system(command)


if __name__ == '__main__':
	main()