import os
import json
import sys
import multiprocessing as mp
import subprocess
import numpy as np

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
	#make new output folders
	# curr_folder = os.getcwd() + '/'


	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	
	folders = []

	job_set_name = "lognorm_radius_test"
	job_set_name = "errorckcsvlognorm"
	job_set_name = "overflowerror"
	job_set_name = "errorckh5lognorm"

	job_set_name = "const_relax"
	job_group = "jobsNovus"

	job_set_name = "constrollingfric"
	job_group = "jobs"
	
	rsize = job_set_name.split("_")[0]

	attempts = [i for i in range(30)]
	# attempts = [23]

	# N = [30,100,300]
	N = [300]
	# M = [3,5,10,15]
	# M=[]
	Temps = [30,100,300,1000]
	# Temps = [30]



	#attempts needs to be first, otherwise follow the order in the list of patterns (called patterns)
	unrolled = u.unroll(attempts,N,Temps)

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 5


	for values in unrolled:
		attempt = values[0]
		n = values[1]
		Temp = values[2]

		job_name = f"relax=1,a={attempt},n={n},t={Temp}"

		#load default input file
		with open(project_path+"default_files/default_input.json",'r') as fp:
			default_input_json = json.load(fp)
		
		# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
		job = default_input_json["data_directory"] + f'{job_group}/' + job_set_name + "relax" + str(attempt) + '/'\
					+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
		copyjob = default_input_json["data_directory"] + f'{job_group}/' + job_set_name + str(attempt) + '/'\
					+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

		

		if os.path.exists(copyjob+"timing.txt") and u.find_max_index(copyjob) == n:
			with open(copyjob+"input.json",'r') as fp:
				input_json = json.load(fp)
			
			if not os.path.exists(job):
				os.makedirs(job)
			else:
				print(f"Job already exists: {job}")


			if os.path.exists(job+"timing.txt"):
				print(f"Sim already complete: {job}")

			elif u.on_queue(job_name):
				print(f"Sim already on queue: {job}")
			else:
				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['N'] = n
				input_json['OMPthreads'] = threadsPerTask
				input_json['MPInodes'] = 1
				input_json['simType'] = "relax"

				# input_json['seed'] = rand_int()
				# input_json['radiiDistribution'] = 'logNormal'
				# input_json['h_min'] = 0.5
				# input_json['dataFormat'] = "csv"
				input_json['relaxIndex'] = n
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

				
				os.system(f"cp {copyjob}{n}* {job}.")
				# os.system(f"cp {copyjob}input.json {job}.")


				sbatchfile = ""
				sbatchfile += "#!/bin/bash\n"
				# sbatchfile += "#SBATCH -C gpu\n"
				sbatchfile += f"#SBATCH -J {job_name}\n"
				sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
				sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
				sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
				sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
				sbatchfile += 'module load hdf5/1.10.8\n'
				
				# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"



				
				with open(job+"sbatch.bash",'w') as sfp:
					sfp.write(sbatchfile)


				#add run script and executable to folders
				# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
				os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
				os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

				folders.append(job)

		elif os.path.exists(copyjob+"timing.txt") and not u.find_max_index(copyjob) == n:
			print(f"!!!!!!!!!!!!!!origin job 'finished' but does not have n={n} file: {copyjob}")
		else:
			print(f"origin job not finished or doesn't exist: {copyjob}")



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

	
