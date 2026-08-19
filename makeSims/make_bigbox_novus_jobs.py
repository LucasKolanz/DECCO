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
		
	job_set_name = "JKRTestest"
	job_set_name = "bigbox"
	# folder_name_scheme = "T_"

	runs_at_once = 1
 
	# attempts = [0] 
	attempts = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14] 
	N = [1000]
	Temps = [3]
	# box_sizes = [1e-3,3e-3,7e-3,1e-2,3e-2]
	box_sizes = [3e-3]
	folders = []

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 28

	for attempt in attempts:
		for n in N:
			for b in box_sizes:
				for Temp in Temps:
					job_name = f"bigboxa={attempt},n={n},b={b:.3e},t={Temp}"
					#load default input file
					with open(project_path+"default_files/default_input.json",'r') as fp:
						input_json = json.load(fp)
					
					# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
					# job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
					# job = input_json["data_directory"] + 'jobs/' + job_set_name + '_' +str(attempt) + '/'
					job = input_json["data_directory"] + f'jobs/{job_set_name}_{attempt}/N_{n}/B_{b:.3e}/T_{Temp}/'
					
					if not os.path.exists(job):
						os.makedirs(job)
					else:
						print(f"Job already exists: {job}")


					if os.path.exists(job+"timing.txt"):
						print(f"Sim already complete: {job}")

					elif u.on_queue(job_name):
						print(f"Sim already on queue: {job}")
					else:
						print(f"(Re)Starting job: {job}")
						####################################
						######Change input values here######
						input_json['temp'] = Temp
						input_json['dynamicTime'] = True
						# input_json['dynamicTime'] = False
						# input_json['seed'] = u.rand_int()
						input_json['seed'] = 101
						input_json['radiiDistribution'] = 'constant'
						# input_json['radiiDistribution'] = 'lognormal'
						input_json['simType'] = 'bigbox'
						# input_json['simType'] = 'custom'
						# input_json['boxDims'] = '1e-4,1e-4,1e-4'
						# input_json['boxDims'] = '5e-4,5e-4,5e-4'
						input_json['boxDims'] = f'{b},{b},{b}'
						# input_json['boxDims'] = '1.0,1.0,1.0'
						# input_json['simType'] = 'custom'
						input_json['N'] = n
						input_json['output_folder'] = job
						# input_json['dataFormat'] = "h5"
						input_json['dataFormat'] = "csv"
						input_json['cor'] = (0.4)**(0.5)
						input_json['OMPthreads'] = threadsPerTask
						input_json['MPInodes'] = 1
						# input_json['impactParameter'] = -1.0
						# input_json['simTimeSeconds'] = 1e-4
						# input_json['timeResolution'] = 1e-6
						# input_json['simTimeSeconds'] = 10e-5
						# input_json['simTimeSeconds'] = 1e-4
						# input_json['simTimeSeconds'] = 3e-3
						input_json['simTimeSeconds'] = 1e-5
						input_json['timeResolution'] = 1e-3
						# input_json['timeResolution'] = 1e-5
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

						#add run script and executable to folders
						# os.system("cp default_files/run_sim.py {}run_sim.py".format(job))
						os.system("cp Collider/Collider.x {}Collider.x".format(job))
						folders.append(job)


						sbatchfile = ""
						sbatchfile += "#!/bin/bash\n"
						# sbatchfile += "#SBATCH -A m2651\n"
						# sbatchfile += "#SBATCH -C gpu\n"
						# sbatchfile += "#SBATCH -q regular\n"
						# sbatchfile += "#SBATCH -t 0:10:00\n"
						sbatchfile += f'#SBATCH --account=lazzati\n'
						sbatchfile += f'#SBATCH --partition=lazzati.q\n'

						#NAME ORDER needs to be same as the file path order
						sbatchfile += f"#SBATCH -J {job_name}\n"
						sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
						sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
						sbatchfile += f"#SBATCH --threads-per-core 1\n"
						sbatchfile += f"#SBATCH --hint=nomultithread\n"
						sbatchfile += f"#SBATCH --exclusive\n"
						sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
						# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

						# sbatchfile += "#SBATCH -G {}\n".format(node)
						# sbatchfile += 'module load gpu\n'

						sbatchfile += f'export OMP_NUM_THREADS={threadsPerTask}\n'
						sbatchfile += f'export OMP_PLACES=cores\n'
						sbatchfile += f'export OMP_PROC_BIND=close\n'
						# sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
						# sbatchfile += 'module load hdf5/1.14.3\n'
						sbatchfile += 'module load gnu12/12.3.0\n'
						sbatchfile += 'module load hdf5/1.10.8\n'
						sbatchfile += 'module load openmpi4/4.1.6\n'
						# sbatchfile += 'module swap openmpi4/4.1.6 mpich\n'

						
						# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"srun --ntasks-per-node={MPITasksPerNode} --cpus-per-task={threadsPerTask} --cpu-bind=socket numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"mpirun --bind-to socket --map-by node numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"mpirun -n {totalMPITasks} numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"


						
						with open(job+"sbatchMulti.bash",'w') as sfp:
							sfp.write(sbatchfile)

		


	print(folders)

	# if args.run:
		# with mp.Pool(processes=runs_at_once) as pool:
		# 	for folder in folders:
		# 		# input_data = inputs[i:i+runs_at_once]
		# 		pool.apply_async(run_job, (folder,))

		# 	pool.close()
		# 	pool.join()

	if args.run:
		cwd = os.getcwd()
		for folder in folders:
			os.chdir(folder)
			os.system('sbatch sbatchMulti.bash')
		os.chdir(cwd)


	
