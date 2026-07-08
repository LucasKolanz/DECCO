import os
import json
import multiprocessing as mp
import subprocess
import argparse
import os
import sys

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u

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

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "bigboxtest"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	# attempts = [0,1,2,3,4,5] 
	attempts = [3] 
	N = [1000]
	Temps = [1000]
	box_sizes = [[]]
	folders = []
	totalMPITasks = 1
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				#load default input file
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)
				
				# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
				# job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				# job = input_json["data_directory"] + 'jobs/' + job_set_name + '_' +str(attempt) + '/'
				job = input_json["data_directory"] + f'jobs/{job_set_name}_{attempt}/N_{n}/T_{Temp}/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))

				# os.system("cp {}/jobs/collidable_aggregate/* {}".format(curr_folder,job))

				####################################
				######Change input values here######
				input_json['temp'] = Temp
				# input_json['seed'] = u.rand_int()
				input_json['seed'] = 101
				input_json['radiiDistribution'] = 'constant'#'lognormal'
				input_json['simType'] = 'bigbox'
				# input_json['boxDims'] = '1e-4,1e-4,1e-4'
				# input_json['boxDims'] = '5e-4,5e-4,5e-4'
				input_json['boxDims'] = '15e-4,15e-4,15e-4'
				# input_json['simType'] = 'custom'
				input_json['N'] = n
				input_json['output_folder'] = job
				# input_json['dataFormat'] = "h5"
				input_json['dataFormat'] = "csv"
				input_json['cor'] = (0.4)**(0.5)
				input_json["OMPthreads"] = 1
				# input_json['impactParameter'] = -1.0
				# input_json['simTimeSeconds'] = 1e-4
				# input_json['timeResolution'] = 1e-6
				# input_json['simTimeSeconds'] = 10e-5
				# input_json['simTimeSeconds'] = 2e-5
				input_json['simTimeSeconds'] = 3e-3
				# input_json['simTimeSeconds'] = 1e-5
				input_json['timeResolution'] = 1e-5
				# input_json['material'] = "amorphousCarbon"
				# input_json['material'] = "quartz"
				# input_json['JKR'] = "False"
				# input_json['JKR'] = "True"
				input_json['density'] = 2.6
				# input_json['relaxIndex'] = n
				# input_json['h_min'] = 0.1
				input_json['h_min'] = 0.5
				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				# input_json['projectileName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				# input_json['targetName'] = "299_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
				# input_json['note'] = "testing"
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#NOVUS CLUSTER
				sbatchfile = ""
				sbatchfile += "#!/bin/bash\n"
				# sbatchfile += "#SBATCH -C gpu\n"
				# sbatchfile += "#SBATCH -q regular\n"
				# sbatchfile += "#SBATCH -t 0:10:00\n"
				# sbatchfile += f'#SBATCH --partition=dri.q\n'

				#FOR ENGR CLUSTER
				#NAME ORDER needs to be same as the file path order
				sbatchfile += f"#SBATCH -J LBGPU,GPUS={totalMPITasks},a={attempt},n={n},t={Temp}\n"
				sbatchfile += "#SBATCH --partition=preempt-gpu.q\n"
				sbatchfile += "#SBATCH --gres=gpu:a40:1\n"
				sbatchfile += f"#SBATCH --ntasks=1\n"
				# sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
				# sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
				# sbatchfile += f"#SBATCH --gpus-per-node=1\n"
				sbatchfile += f"#SBATCH --cpus-per-task 1\n\n"
				# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

				# sbatchfile += "#SBATCH -G {}\n".format(node)
				# sbatchfile += 'module load gpu\n'

				# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
				sbatchfile += 'lscpu\n'
				# sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
				# sbatchfile += 'module load hdf5/1.14.3\n'
				sbatchfile += 'module load hdf5/1.10.8\n'
				# sbatchfile += 'module load gcc/8.3\n'
				# sbatchfile += 'module load gnu12/12.3.0\n'

				# sbatchfile += 'module load mpich/3.3\n'
				sbatchfile += 'module list\n'
				sbatchfile += 'echo "Running on host: $(hostname)"\n'
				# sbatchfile += 'export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH\n'
				# sbatchfile += 'echo $LD_LIBRARY_PATH\n'
				# sbatchfile += 'ls /usr/lib64/libc.*\n'
				# sbatchfile += "unset NV_ACC_DEBUG\n"
				# sbatchfile += "unset NV_ACC_NOTIFY\n"
				# sbatchfile += "export NV_ACC_NOTIFY=3\n"
				# sbatchfile += "export NV_ACC_DEBUG=0x800\n"
				sbatchfile += "export ACC_DEVICE_TYPE=nvidia\n"
				sbatchfile += "echo $LD_LIBRARY_PATH\n"
				sbatchfile += "export HDF5_USE_FILE_LOCKING=FALSE\n"
				sbatchfile += "ls /opt/ohpc/pub/libs/gnu12/hdf5/1.10.8/lib/\n"
				sbatchfile += "nvidia-smi\n"

			

				
				# sbatchfile += f"srun -n {totalNodes} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				# sbatchfile += f"srun --ntasks-per-node={MPITasksPerNode} --cpus-per-task={threadsPerTask} --cpu-bind=socket numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				# sbatchfile += f"mpirun --bind-to socket --map-by node numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				sbatchfile += f"mpirun -np 1 {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
				# sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"


				
				with open(job+"sbatchMulti.bash",'w') as sfp:
					sfp.write(sbatchfile)



				#add run script and executable to folders
				# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
				os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
				os.system(f"cp {project_path}Collider/ball_group.cpp {job}ball_group.cpp")
				os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

				
					
				folders.append(job)
	# print(folders)


	print(folders)
	if args.run:
		cwd = os.getcwd()
		for folder in folders:
			os.chdir(folder)
			os.system('sbatch sbatchMulti.bash')
		os.chdir(cwd)

	