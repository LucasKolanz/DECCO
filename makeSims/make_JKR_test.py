import os
import json
import multiprocessing as mp
import subprocess
import os

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


# def rand_int():
# 	# Generating a random integer from 0 to the maximum unsigned integer in C++
# 	# In C++, the maximum value for an unsigned int is typically 2^32 - 1
# 	max_unsigned_int_cpp = 2**32 - 1
# 	random_unsigned_int = random.randint(0, max_unsigned_int_cpp)
# 	return random_unsigned_int

def run_job(location):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = [f"{location}Collider.x",location]

	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)

if __name__ == '__main__':
	#make new output folders
	# curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "JKRTest"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [0] 
	N = [10]
	Temps = [1]
	folders = []
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				#load default input file
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)
				
				# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
				# job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				job = input_json["data_directory"] + 'jobs/' + job_set_name + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))

				# os.system("cp {}/jobs/collidable_aggregate/* {}".format(curr_folder,job))

				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['seed'] = 101
				input_json['radiiDistribution'] = 'constant'#'lognormal'
				input_json['simType'] = 'Collider'
				input_json['N'] = n
				input_json['output_folder'] = job
				input_json['dataFormat'] = "csv"
				input_json['impactParameter'] = -1.0
				input_json['simTimeSeconds'] = 1e-4
				input_json['timeResolution'] = 1e-5
				input_json['material'] = "amorphousCarbon"
				input_json['JKR'] = "True"
				# input_json['kConsts'] = 3e3
				# input_json['relaxIndex'] = n-3
				input_json['h_min'] = 0.5
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
				folders.append(job)
	


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

	# print(folders)
	# cwd = os.getcwd()
	# for folder in folders:
	# 	os.chdir(folder)
	# 	os.system('qsub qsub.bash')
	# os.chdir(cwd)


	
