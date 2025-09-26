import os
import json
import multiprocessing as mp
import subprocess
import numpy as np

def run_job(location,num_balls):
	cmd = ["python3", "{}run_multicore_sim.py".format(location), location, str(num_balls)]
	# print(cmd)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","ColliderMultiCore"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	# job_set_name = "openMPallLoops"
	# job_set_name = "strongScaleGrow"
	# job_set_name = "weakScaleCollide"
	job_set_name = "weakScaleCollide_O2_"
	# job_set_name = "pipeAndOpenmp"
	# job_set_name = "smallerDt"
	# job_set_name = "forceTest"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [2,3,4,5]
	attempts = [2]
	threads = [1,4,16,64]
	genBalls = [76,151,301,601]
	N = [1]
	Temps = [100]
	folders = []
	for attempt in attempts:
		for n in N:
			for t_i,thread in enumerate(threads):
				# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
				job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'thread_' + str(thread) + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))

				# os.system("cp {}/jobs/collidable_aggregate/* {}".format(curr_folder,job))

				#load default input file
				with open(curr_folder+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				####################################
				######Change input values here######
				input_json['temp'] = Temps[0]
				input_json['seed'] = 101
				input_json['radiiDistribution'] = 'constant'
				# input_json['kConsts'] = 3e3
				input_json['N'] = n
				input_json['simType'] = "BPCA"
				input_json['genBalls'] = genBalls[t_i]
				input_json['h_min'] = 0.5
				input_json['OMPthreads'] = thread
				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				input_json['projectileName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				input_json['targetName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				# input_json['note'] = "Uses openmp and loop unwinding to parallelize sim_one_step."
				input_json['note'] = "Weak scaling with {} OMP threads.".format(thread)
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)


				# time = np.ceil(n/15)

				sbatchfile = ""
				sbatchfile += "#!/bin/bash\n"
				sbatchfile += "#SBATCH -A m2651\n"
				sbatchfile += "#SBATCH -C cpu\n"
				sbatchfile += "#SBATCH -q regular\n"
				sbatchfile += "#SBATCH -t 3:00:00\n"
				sbatchfile += "#SBATCH -N 1\n"
				# sbatchfile += "#SBATCH -c {}\n\n".foramt(2*thread)
				# sbatchfile += 'module load gpu\n'
				# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(thread)
				sbatchfile += 'export SLURM_CPU_BIND="cores"\n'
				sbatchfile += "srun -n 1 -c {} --cpu-bind=cores ./ColliderMultiCore.x {} 2>sim_err.log 1>sim_out.log".format(thread*2,job)
				
				with open(job+"sbatchMulti.bash",'w') as sfp:
					sfp.write(sbatchfile)

				#add run script and executable to folders
				os.system("cp default_files/run_multicore_sim.py {}run_multicore_sim.py".format(job))
				os.system("cp ColliderMultiCore/ColliderMultiCore.x {}ColliderMultiCore.x".format(job))
				os.system("cp ColliderMultiCore/ColliderMultiCore.cpp {}ColliderMultiCore.cpp".format(job))
				os.system("cp ColliderMultiCore/ball_group_multi_core.hpp {}ball_group_multi_core.hpp".format(job))
				# os.system("cp jobs/collidable_aggregate/* {}".format(job))

				folders.append(job)
	# print(folders)
	# if len(N) != len(folders):
	# 	N = [str(N[0]) for i in range(len(folders))]

	# inputs = list(zip(folders,N))
	print(folders)
	
	cwd = os.getcwd()
	print(folders)
	for folder in folders:
		os.chdir(folder)
		os.system("sbatch sbatchMulti.bash")
	os.chdir(cwd)

	# for i in range(0,len(folders),runs_at_once):
	# 	with mp.Pool(processes=runs_at_once) as pool:
	# 		pool.starmap(run_job,folders[i:i+runs_at_once]) 


	
