import os
import json
import multiprocessing as mp
import subprocess
import random
from datetime import datetime
import re


relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


random.seed(int(datetime.now().timestamp()))

def get_squeue_output():
    try:
        # Run the squeue command and capture its output

        result = subprocess.run(['squeue', '-o', '"%.20u %.25j"'], capture_output=True, text=True)
        output = result.stdout
        return output
    except subprocess.CalledProcessError as e:
        # Handle any errors that occur during the command execution
        print(f"Error executing squeue: {e}")
        return None


def run_job(location):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = [f"{location}Collider.x",location]

	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)


def same_job(fullpath, job_name):

	fpsplit = fullpath.split('/')

	fpattrs = re.split(r'\D+',"".join(fpsplit[-4:-1]))
	fpattrs = [int(i) for i in fpattrs if len(i) > 0]
	
	qattrs = re.split(r'\D+',job_name)
	qattrs = [int(i) for i in qattrs if len(i) > 0]

	if len(fpattrs) != len(qattrs):
		print("ERROR IN same_job")
		exit(0)

	for i in range(len(qattrs)):
		if fpattrs[i] != qattrs[i]:
			return False
	return True

def on_queue(fullpath):
	queue_out = get_squeue_output()
	for line in queue_out.split('\n')[1:]:
		line = line.strip('"').split()
		if len(line) > 0:
			if line[0] == "kolanzl":
				if same_job(fullpath,line[1]):
					return True
	return False

def rand_int(start=0,stop=-1):
	# Generating a random integer from 0 to the maximum unsigned integer in C++
	# In C++, the maximum value for an unsigned int is typically 2^32 - 1
	if stop < 0:
		max_unsigned_int_cpp = 2**32 - 1
		stop = max_unsigned_int_cpp
	random_unsigned_int = random.randint(start, stop)
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
	job_set_name = "lognorm_collisions"
	job_set_name = "const_collisions"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	

	# attempts = [i for i in range(30)]
	attempts = [i for i in range(10)]
	attempts = [0]
	# attempts = [90]

	# attempts_300 = attempts

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 1
	# N_combos = [[30,30],[30,100],[30,300]]
	N_combos = [[30,100]]
	# N = [300]
	etas = [i/2.0+0.5 for i in range(0,24)]
	etas = [etas[0]]


	# attempts_300 = attempts
	# node = 1
	folders = []
	for ns in N_combos:
		# threads = 1
		# if n == 30:
		# 	threads = 1
		# elif n == 100:
		# 	threads = 2
		# else:# n == 300:
		# 	threads = 24
		# temp_attempt = attempts
		# if n == 300:
		# 	temp_attempt = attempts_300
		for attempt in attempts:
			for eta in etas:
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(ns[0])+'-'+str(ns[1]) + '/' + 'eta_' + str(eta) + '/'

				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				if os.path.exists(job+"timing.txt"):
					print("Sim already complete")

				# elif on_queue(job):
				# 	print("Sim already on queue")
				else:
					#load default input file

					####################################
					######Change input values here######
					# input_json['temp'] = Temp
					input_json['seed'] = rand_int()
					# input_json['radiiDistribution'] = 'logNormal'
					input_json['simType'] = 'collider'
					input_json['eta'] = eta
					# input_json['N'] = n
					input_json['h_min'] = 0.5
					input_json['dataFormat'] = "h5"
					input_json['output_folder'] = job
					input_json['OMPthreads'] = threadsPerTask
					input_json['MPInodes'] = totalMPITasks
					input_json['v_custom'] = 1
					
					input_json['projectileName'] = f"/media/kolanzl/easystore/SpaceLab_data/jobsNovus/{job_set_name.split('_')[0]}_relax{rand_int(0,29)}/N_{ns[0]}/T_{'3'}/{ns[0]-3}_RELAXdata.h5"
					input_json['targetName'] = f"/media/kolanzl/easystore/SpaceLab_data/jobsNovus/{job_set_name.split('_')[0]}_relax{rand_int(0,29)}/N_{ns[1]}/T_{3}/{ns[1]-3}_RELAXdata.h5"
		
					# input_json['u_s'] = 0.5
					# input_json['u_r'] = 0.5
					input_json['note'] = "For colliding BPCA grown aggregates."
					####################################

					with open(job + "input.json",'w') as fp:
						json.dump(input_json,fp,indent=4)

					# sbatchfile = ""
					# sbatchfile += "#!/bin/bash\n"
					# # sbatchfile += "#SBATCH -A m2651\n"
					# # sbatchfile += "#SBATCH -C gpu\n"
					# # sbatchfile += "#SBATCH -q regular\n"
					# # sbatchfile += "#SBATCH -t 0:10:00\n"
					# # sbatchfile += f'#SBATCH --partition=dri.q\n'
					# sbatchfile += f"#SBATCH -J a={attempt},n={ns[0]}-{ns[1]},e={eta}\n"
					# sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
					# sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
					# sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
					# # sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

					# # sbatchfile += "#SBATCH -G {}\n".format(node)
					# # sbatchfile += 'module load gpu\n'

					# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
					# # sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
					# # sbatchfile += 'module load hdf5/1.14.3\n'
					# # sbatchfile += 'module load hdf5/1.10.8\n'
					# sbatchfile += 'module load gnu12/12.3.0\n'
					# sbatchfile += 'module load openmpi4/4.1.6\n'
					# # sbatchfile += 'module swap openmpi4/4.1.6 mpich\n'

					
					# # sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
					# # sbatchfile += f"srun --ntasks-per-node={MPITasksPerNode} --cpus-per-task={threadsPerTask} --cpu-bind=socket numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
					# # sbatchfile += f"mpirun --bind-to socket --map-by node numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
					# # sbatchfile += f"mpirun -n {totalMPITasks} numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
					# sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"

					# with open(job+"sbatchMulti.bash",'w') as sfp:
					# 	sfp.write(sbatchfile)

					

					#add run script and executable to folders
					# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
					os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
					os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
					os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

					folders.append(job)

	print(folders)
	# cwd = os.getcwd()
	# for folder in folders:
	# 	os.chdir(folder)
	# 	os.system('sbatch sbatchMulti.bash')
	# os.chdir(cwd)

	with mp.Pool(processes=runs_at_once) as pool:
		for folder in folders:
			# input_data = inputs[i:i+runs_at_once]
			pool.apply_async(run_job, (folder,))

		pool.close()
		pool.join()








	
