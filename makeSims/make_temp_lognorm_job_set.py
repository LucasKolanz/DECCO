import os
import json
import multiprocessing as mp
import subprocess
import random
import argparse
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

	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)


	job_set_name = "overflow_tester"
	job_set_name = "test_"
	job_set_name = "lognorm_"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	# attempts = [i for i in range(10)]
	attempts = [7,19,20] 
	# attempts_300 = [i for i in range(5)]

	#test it out first
	# attempts = [0]
	# attempts_300 = [0]

	# N = [30,100,300]
	N = [300]
	# N = [5]
	# Temps = [3,10,30,100,300,1000]
	Temps = [1000]

	folders = []

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_template = input_json["data_directory"] + 'jobsCosine/' + job_set_name + '{a}/N_{n}/T_{t}/'


	for n in N:
		for Temp in Temps:
			# temp_attempt = attempts
			# if n == 300:
			# 	temp_attempt = attempts_300
			# for attempt in temp_attempt:
			for attempt in attempts:
				job = job_template.replace('{a}',str(attempt)).replace('{n}',str(n)).replace('{t}',str(Temp))
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))

				if os.path.exists(job+'timing.txt'):
					print("Job already finished.")
				else:
					#load default input file
					with open(project_path+"default_files/default_input.json",'r') as fp:
						input_json = json.load(fp)

					####################################
					######Change input values here######
					input_json['JKR'] = "false"
					input_json['temp'] = Temp
					input_json['N'] = n
					input_json['seed'] = u.rand_int()
					input_json['radiiDistribution'] = 'logNormal'
					input_json['impactParameter'] = -1.0
					input_json['OMPthreads'] = 2
					input_json['h_min'] = 0.5
					input_json['dataFormat'] = "csv"
					# input_json['dataFormat'] = "h5"
					input_json['output_folder'] = job
					input_json['note'] = "Finally running these"
					####################################

					with open(job + "input.json",'w') as fp:
						json.dump(input_json,fp,indent=4)


					
					#####################3#add run script and executable to folders
					os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
					# os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
					# os.system(f"cp {project_path}Collider/ball_group.cpp {job}ball_group.cpp")
					# os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

					folders.append(job)
					######################################################


	print(folders)

	if args.run:
		with mp.Pool(processes=runs_at_once) as pool:
			for folder in folders:
				pool.apply_async(run_job,(folder,)) 

			pool.close()
			pool.join()




	
