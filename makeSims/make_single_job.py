import os
import json
import multiprocessing as mp
import subprocess

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


	# out = os.system("./ColliderSingleCore.o {}".format(curr_folder))
	# out = os.system("./ColliderSingleCore.o {} 1>> {} 2>> {}".format(curr_folder,output_file,error_file))
	
	# cmd = ["srun","-n","1","-c","2","{}ColliderSingleCore.x".format(location), location, str(num_balls)]

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
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "lognorm_radius_test"
	job_set_name = "errorckcsvlognorm"
	job_set_name = "errorckh5lognorm"
	job_set_name = "overflowerror"
	job_set_name = "TESTBCCA"
	# job_set_name = "TEST"

	# folder_name_scheme = "T_"

	SPECIAL_FOLDER = ""#"/home/lucas/Desktop/SpaceLab_data/largejob/"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [0] 
	N = [8]
	threads = []
	# Temps = [3,10,30,100,300,1000]
	Temps = [3]
	folders = []
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				#load default input file
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)
				
				# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
				job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))

				####################################
				######Change input values here######
				# input_json['temp'] = Temp
				input_json['N'] = n
				input_json['output_folder'] = job
				input_json['OMPthreads'] = 1
				input_json['MPInodes'] = 1

				# input_json['seed'] = rand_int()
				input_json['seed'] = 100
				input_json['radiiDistribution'] = 'logNormal'
				# input_json['h_min'] = 0.5
				input_json['dataFormat'] = "csv"
				input_json['simType'] = "BCCA"
				input_json['symmetric'] = "true"
				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				# input_json['note'] = "Does this work at all?"
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
				os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
				os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

				# os.system(f"cp /home/lucas/Desktop/SpaceLab_data/test2/N_5/T_3/*data.h5 {job}.")
				
				folders.append(job)
	# print(folders)


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

	
