import os
import json
import multiprocessing as mp
import subprocess
import random
from datetime import datetime


relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

random.seed(datetime.now().timestamp())

def rand_int():
	# Generating a random integer from 0 to the maximum unsigned integer in C++
	# In C++, the maximum value for an unsigned int is typically 2^32 - 1
	max_unsigned_int_cpp = 2**32 - 1
	random_unsigned_int = random.randint(0, max_unsigned_int_cpp)
	return random_unsigned_int

def run_job(location):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = [f"{location}Collider.x",location]

	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	print(project_path)
	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "GPUtest"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [1] 
	N = [1203]
	Temps = [10]
	folders = []
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				#load default input file
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)
				job = input_json['data_directory'] + 'jobs/' + job_set_name + str(attempt) + '/'
				# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
							# + 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))
				# os.system("cp {}/jobs/collidable_aggregate/* {}".format(curr_folder,job))


				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['seed'] = 101
				input_json['radiiDistribution'] = 'constant'
				# input_json['kConsts'] = 3e3
				input_json['h_min'] = 0.5
				input_json['dataFormat'] = "csv"
				input_json['N'] = n
				input_json['output_folder'] = job
				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				# input_json['projectileName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				# input_json['targetName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				input_json['note'] = "testing"
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				# os.system("cp default_files/run_sim.py {}run_sim.py".format(job))
				os.system("cp Collider/Collider.x {}Collider.x".format(job))
				os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
				os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

				os.system(f"cp /mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_data/jobsCosine/lognorm11/N_300/T_10/10_* {job}.")
				os.system(f"touch {job}11_constants.csv")
				os.system(f"touch {job}11_simData.csv")
				os.system(f"touch {job}11_energy.csv")
				folders.append(job)
	

	print(folders)

	with mp.Pool(processes=runs_at_once) as pool:
		for folder in folders:
			# input_data = inputs[i:i+runs_at_once]
			pool.apply_async(run_job, (folder,))

		pool.close()
		pool.join()


	
