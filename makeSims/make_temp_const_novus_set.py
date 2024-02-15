import os
import json
import multiprocessing as mp
import subprocess
import random

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

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

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)


	job_set_name = "TEST"
	job_set_name = "const"
	# folder_name_scheme = "T_"

	# runs_at_once = 7
	# attempts = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	# attempts = [i for i in range(30)]
	attempts = [i for i in range(30)]
	# attempts = [30]
	# attempts_300 = [i for i in range(30)]
	# attempts = [1] 
	attempts_300 = attempts

	#test it out first
	# attempts = [0]
	# attempts_300 = [0]

	node = 1
	N = [30,100,300]
	# N = [300]
	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	folders = []
	for n in N:
		threads = 1
		if n == 30:
			threads = 1
		elif n == 100:
			threads = 2
		else:# n == 300:
			threads = 16
		temp_attempt = attempts
		if n == 300:
			temp_attempt = attempts_300
		for attempt in temp_attempt:
			for Temp in Temps:
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				if os.path.exists(job+"timing.txt"):
					print("Sim already complete")
				else:
					#load default input file

					####################################
					######Change input values here######
					input_json['temp'] = Temp
					input_json['seed'] = rand_int()
					input_json['radiiDistribution'] = 'constant'
					input_json['N'] = n
					input_json['h_min'] = 0.5
					input_json['dataFormat'] = "h5"
					input_json['output_folder'] = job
					input_json['OMPthreads'] = threads
					# input_json['u_s'] = 0.5
					# input_json['u_r'] = 0.5
					input_json['note'] = "Rerunning constant size ball runs."
					####################################

					with open(job + "input.json",'w') as fp:
						json.dump(input_json,fp,indent=4)


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


					
					with open(job+"sbatchMulti.bash",'w') as sfp:
						sfp.write(sbatchfile)

					#add run script and executable to folders
					# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
					os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
					os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
					os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

					folders.append(job)

print(folders)
cwd = os.getcwd()
for folder in folders:
	os.chdir(folder)
	os.system('sbatch sbatchMulti.bash')
os.chdir(cwd)







	
