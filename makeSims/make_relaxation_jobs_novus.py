import os
import json
import multiprocessing as mp
import subprocess
import numpy as np

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
	job_set_name = "overflowerror"
	job_set_name = "errorckh5lognorm"

	job_set_name = "const_relax"
	job_group = "jobsNovus"

	job_set_name = "lognorm_relax"
	job_group = "jobsCosine"
	job_group = "jobs"
	
	rsize = job_set_name.split("_")[0]

	# folder_name_scheme = "T_"


	# # runs_at_once = 12
	# attempts = [i for i in range(30)] 
	# # attempts = [18] 
	# N = [30,100,300]
	# # N = [30]
	# Temps = [3,10,30,100,300,1000]
	# # Temps = [1000]
	# folders = []
	# # threads = []
	# # for attempt in attempts:
	# # 	for n in N:
	# # 		for Temp in Temps:

	runvals = [[9,300,3],[26,100,1000],[22,300,3],[23,300,3],[24,300,3],[26,300,3],[29,300,3],[22,300,10],[23,300,10],[24,300,10],[25,300,10],[26,300,10],[24,300,30],[25,300,30],[4,300,1000]]
	runvals = np.array(runvals)
	# print(runvals.shape)
	# exit(0)
	for i in range(runvals.shape[0]):
		attempt = runvals[i,0]
		n = runvals[i,1]
		Temp = runvals[i,2]


		#load default input file
		with open(project_path+"default_files/default_input.json",'r') as fp:
			default_input_json = json.load(fp)
		
		# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
		job = default_input_json["data_directory"] + f'{job_group}/' + job_set_name + str(attempt) + '/'\
					+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
		copyjob = default_input_json["data_directory"] + f'{job_group}/' + rsize + str(attempt) + '/'\
					+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
		

		if os.path.exists(copyjob+"timing.txt"):
			if not os.path.exists(job):
				os.makedirs(job)
			else:
				print("Job '{}' already exists.".format(job))

			if not os.path.exists(job+"timing.txt"):
				with open(copyjob+"input.json",'r') as fp:
					input_json = json.load(fp)
				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['N'] = n
				input_json['OMPthreads'] = 1
				input_json['MPInodes'] = 1
				input_json['simType'] = "relax"

				# input_json['seed'] = rand_int()
				# input_json['radiiDistribution'] = 'logNormal'
				# input_json['h_min'] = 0.5
				# input_json['dataFormat'] = "csv"
				input_json['relaxIndex'] = n-3
				input_json['simTimeSeconds'] = 1e-2
				# input_json['timeResolution'] = 

				input_json['output_folder'] = job
				input_json['data_directory'] = default_input_json['data_directory']
				input_json['project_directory'] = default_input_json['project_directory']

				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				input_json['note'] = "Relaxation job"

				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				
				os.system(f"cp {copyjob}{n-3}* {job}.")
				# os.system(f"cp {copyjob}input.json {job}.")


				sbatchfile = ""
				sbatchfile += "#!/bin/bash\n"
				# sbatchfile += "#SBATCH -A m2651\n"
				# sbatchfile += "#SBATCH -C gpu\n"
				# sbatchfile += "#SBATCH -q regular\n"
				# sbatchfile += "#SBATCH -t 0:10:00\n"
				sbatchfile += f"#SBATCH -J {job_set_name}\n"
				sbatchfile += f"#SBATCH -N {node}\n"
				sbatchfile += f"#SBATCH -n {node}\n"
				sbatchfile += f"#SBATCH -c {threads}\n\n"
				# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

				# sbatchfile += "#SBATCH -G {}\n".format(node)
				# sbatchfile += 'module load gpu\n'

				sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threads)
				sbatchfile += 'export SLURM_CPU_BIND="cores"\n'
				# sbatchfile += 'module load hdf5/1.14.3\n'
				sbatchfile += 'module load hdf5/1.10.8\n'
				
				# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"


				
				with open(job+"sbatch.bash",'w') as sfp:
					sfp.write(sbatchfile)


				#add run script and executable to folders
				# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
				os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
				os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

				folders.append(job)
		else:
			print(f"origin job doesn't exist: {copyjob}")

		# print(folders)



	# for i in range(0,len(folders),runs_at_once):
	# 	with mp.Pool(processes=runs_at_once) as pool:
	# 		pool.starmap(run_job,inputs[i:i+runs_at_once]) 
	
	# print(folders)
	# with mp.Pool(processes=runs_at_once) as pool:
	# 	for folder in folders:
	# 		# input_data = inputs[i:i+runs_at_once]
	# 		pool.apply_async(run_job, (folder,))

	# 	pool.close()
	# 	pool.join()

	print(folders)
	# cwd = os.getcwd()
	# for folder in folders:
	# 	os.chdir(folder)
	# 	os.system('sbatch sbatch.bash')
	# os.chdir(cwd)

	
