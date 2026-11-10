import os
import json
import multiprocessing as mp
import subprocess
import re
import argparse
import random
import sys


relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u


def run_job(location):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = [f"{location}Collider.x",location]

	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)

if __name__ == '__main__':

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	parser = argparse.ArgumentParser(
	description="Prepare DEM jobs and optionally submit them via Slurm."
	)
	parser.add_argument(
		"-r",
		"--run",
		action="store_true",
		help="Actually submit jobs with sbatch (otherwise do a dry run).",
	)

	args = parser.parse_args()

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "weldVerificationCollision"

	# folder_name_scheme = "T_"

	# SPECIAL_FOLDER = ""#"/home/lucas/Desktop/SpaceLab_data/largejob/"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [6] 
	folders = []

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 16
	
	for attempt in attempts:
		
		# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
		job = input_json["data_directory"] + 'jobs/' + job_set_name + '_' + str(attempt) + '/'
		
		job_name = f"WELDVerification_a={attempt}"

		if not os.path.exists(job):
			os.makedirs(job)
		else:
			print(f"Job already exists: {job}")


		if os.path.exists(job+"timing.txt"):
			print(f"Sim already complete: {job}")

		print(f"(Re)Starting job: {job}")

		####################################
		######Change input values here######
		# input_json['temp'] = 3
		# input_json['N'] = 5
		input_json['seed'] = 101
		# input_json['weld'] = 'false'
		input_json['weld'] = 'true'
		input_json['JKR'] = "false"
		input_json['output_folder'] = job
		input_json['OMPthreads'] = threadsPerTask
		input_json['MPInodes'] = totalMPITasks

		# input_json['simTimeSeconds'] = 0.001  
		# input_json['timeResolution'] = 1e-5

		input_json['simTimeSeconds'] = 0.001
		input_json['timeResolution'] = 0.001/(1000.0)

		input_json['v_custom'] = '10'
		# input_json['h_min'] = 0.5
		input_json['dataFormat'] = "csv"
		# input_json['simType'] = "BPCA"
		input_json['simType'] = "collider"
		input_json['impactParameter'] = 0.5
		
		input_json['projectileName'] = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_150/N_300/T_1000/150_simData.csv"
		input_json['targetName'] = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_150/N_300/T_1000/150_simData.csv"

		####################################
		with open(job + "input.json",'w') as fp:
			json.dump(input_json,fp,indent=4)

		#add run script and executable to folders 
		# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
		os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
		# os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
		# os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

		folders.append(job)



	###############For running locally
	print(folders)

	# for i in range(0,len(folders),runs_at_once):
	# 	with mp.Pool(processes=runs_at_once) as pool:
	# 		pool.starmap(run_job,inputs[i:i+runs_at_once]) 

	with mp.Pool(processes=runs_at_once) as pool:
		for folder in folders:
			# input_data = inputs[i:i+runs_at_once]
			pool.apply_async(run_job, (folder,))

		pool.close()
		pool.join()

	
