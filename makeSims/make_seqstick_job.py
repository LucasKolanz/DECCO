import os
import json
import multiprocessing as mp
import subprocess
import random

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


	# out = os.system("./ColliderSingleCore.o {}".format(curr_folder))
	# out = os.system("./ColliderSingleCore.o {} 1>> {} 2>> {}".format(curr_folder,output_file,error_file))
	
	# cmd = ["srun","-n","1","-c","2","{}ColliderSingleCore.x".format(location), location, str(num_balls)]

def rand_int():
	# Generating a random integer from 0 to the maximum unsigned integer in C++
	# In C++, the maximum value for an unsigned int is typically 2^32 - 1
	max_unsigned_int_cpp = 2**32 - 1
	random_unsigned_int = random.randint(0, max_unsigned_int_cpp)
	return random_unsigned_int

def run_job(location,mode,n):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = ["python3", f"{location}seqStick.py", "-d", f"{mode}", "-n", f"{n}"]


	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)

if __name__ == '__main__':

		

	# job_set_name = "TESTBAPA"
	job_set_name = "SeqStickLognorm_"

	# folder_name_scheme = "T_"


	runs_at_once = 10
	# attempts = [2] 
	attempts = [i for i in range(0,30)]#[0,1,2,3,4,5,6,7,8,9]#,11,12,13,14,15,16,17,18,19,20] 
	# attempts = [0]
	N = [300] #final size
	folders = []

	mode = "lognormal"

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_template = input_json["data_directory"] + 'jobs/' + job_set_name + '{a}/N_{n}/'

	for attempt in attempts:
		for n in N:
			#load default input file
			# with open(project_path+"default_files/default_input.json",'r') as fp:
			# 	input_json = json.load(fp)
			
			# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
			job = job_template.replace('{a}',str(attempt)).replace('{n}',str(n))

			
			if not os.path.exists(job):
				os.makedirs(job)
			else:
				print(f"Job '{job}' already exists.")

			
			if not os.path.exists(job+f'timing.txt'):
				#add run script and executable to folders
				# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}SequentialSticking/seqStick.py {job}seqStick.py")
				folders.append(job)
			else:
				print(f"Job '{job}' already complete.")



	print(folders)

	
	with mp.Pool(processes=runs_at_once) as pool:
		for folder in folders:
			# input_data = inputs[i:i+runs_at_once]
			pool.apply_async(run_job, (folder,mode,N[0]))

		pool.close()
		pool.join()

	
