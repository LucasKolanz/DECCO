import os
import json
import multiprocessing as mp
import subprocess
import argparse
import random
import glob
import sys

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u

	# out = os.system("./ColliderSingleCore.o {}".format(curr_folder))
	# out = os.system("./ColliderSingleCore.o {} 1>> {} 2>> {}".format(curr_folder,output_file,error_file))
	
	# cmd = ["srun","-n","1","-c","2","{}ColliderSingleCore.x".format(location), location, str(num_balls)]


def run_job(location):
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	cmd = [f"{location}Collider.x",location]

	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Prepare DEM jobs and optionally submit them."
	)
	parser.add_argument(
		"-r",
		"--run",
		action="store_true",
		help="Actually submit jobs with sbatch (otherwise do a dry run).",
	)

	args = parser.parse_args()

	#make new output folders
	# curr_folder = os.getcwd() + '/'


	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		

	job_set_name = "TESTASYMBAPA"
	job_set_name = "AsymBAPA"

	# folder_name_scheme = "T_"


	# runs_at_once = 1
	# attempts = [2] 
	# attempts = [i for i in range(0,10)]
	attempts = [i for i in range(0,5)]
	attempts = [i for i in range(5,10)]
	# attempts = [8]
	N = [300] #minimum final size
	# M = [20,30,50,60,100] #starting sizes
	M = [-1] #starting sizes
	threads = []
	# Temps = [3,10,30,100,300,1000]
	Temps = [1000]
	folders = []

	m_range = [1,100]

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_template = input_json["data_directory"] + 'jobs/' + job_set_name + '_{a}/N_{n}/T_{t}/'

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 32

	for attempt in attempts:
		for m in M:
			for n in N:
				for Temp in Temps:
					#load default input file
					# with open(project_path+"default_files/default_input.json",'r') as fp:
					# 	input_json = json.load(fp)
					
					# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
					job = job_template.replace('{a}',str(attempt)).replace('{n}',str(n)).replace('{t}',str(Temp))

					job_name = f"ASBA=1,a={attempt},m={m},n={n},t={Temp}"
					
					if not os.path.exists(job):
						os.makedirs(job)
					else:
						print("Job '{}' already exists.".format(job))

					if os.path.exists(job+"timing.txt"):
						print(f"Sim already complete: {job}")

					elif u.on_queue(job_name):
						print(f"Sim already on queue: {job}")
					else:
						

						####################################
						######Change input values here######
						input_json['JKR'] = 'false'
						input_json['temp'] = Temp
						input_json['N'] = n
						input_json['M'] = m
						input_json['M_range'] = f"{m_range[0]}:{m_range[1]}"
						input_json['output_folder'] = job
						input_json['OMPthreads'] = threadsPerTask
						input_json['MPInodes'] = 1
						input_json['impactParameter'] = -1.0

						input_json['seed'] = u.rand_int()
						# input_json['seed'] = 101
						input_json['genBalls'] = 2

						input_json['radiiDistribution'] = 'logNormal'
						# input_json['h_min'] = 0.5
						
						input_json['timeResolution'] = 1e-5
						input_json['simTimeSeconds'] = 5e-4

						# input_json['simTimeSeconds'] = 5e-4

						input_json['dataFormat'] = "csv"
						input_json['simType'] = "BAPA"
						input_json['random_folder_template'] = input_json['data_directory']+"localLognormData/lognorm{a}/N_30/T_1000/"

						# input_json['u_s'] = 0.5
						# input_json['u_r'] = 0.5
						# input_json['note'] = "Does this work at all?"
						####################################

						print("dumping json")
						with open(job + "input.json",'w') as fp:
							json.dump(input_json,fp,indent=4)

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
						sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
						# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

						# sbatchfile += "#SBATCH -G {}\n".format(node)
						# sbatchfile += 'module load gpu\n'

						sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
						# sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
						# sbatchfile += 'module load hdf5/1.14.3\n'
						sbatchfile += 'module load gnu12/12.4.0\n'
						sbatchfile += 'module load hdf5/1.10.8\n'
						# sbatchfile += 'module load openmpi4/4.1.6\n'
						# sbatchfile += 'module swap openmpi4/4.1.6 mpich\n'

						
						# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"srun --ntasks-per-node={MPITasksPerNode} --cpus-per-task={threadsPerTask} --cpu-bind=socket numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"mpirun --bind-to socket --map-by node numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						# sbatchfile += f"mpirun -n {totalMPITasks} numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
						sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"

						print("dumping sbatch file")
						with open(job+"sbatchMulti.bash",'w') as sfp:
							sfp.write(sbatchfile)





						print("copying executable")
						#add run script and executable to folders
						# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
						os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
						# os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
						# os.system(f"cp {project_path}Collider/ball_group.cpp {job}ball_group.cpp")
						# os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

						# randm = random.randint(m_range[0],m_range[1])
						randm = 3
						if randm >= 3 and not os.path.exists(f"{job}/{randm}_checkpoint.txt"): 
							print("getting random source")

							for _ in range(100):
								randint = random.randint(0, 29)
								source = f"{input_json['data_directory']}/localLognormData/lognorm{randint}/N_300/T_{Temp}/{randm}_*"
								if any(glob.iglob(source)):
									break
							else:
								raise FileNotFoundError("No matching file found after many attempts.")

							os.system(f"cp {source} {job}")
							os.system(f"touch {source} {job}/{randm}_checkpoint.txt")
								# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm{randint}/N_300/T_3/{m}_constants.csv {job}{m}_constants.csv")
								# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm{randint}/N_300/T_3/{m}_simData.csv {job}{m}_simData.csv")
								# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm{randint}/N_300/T_3/{m}_energy.csv {job}{m}_energy.csv")
							
						folders.append(job)



	print(folders)

	if args.run:
		cwd = os.getcwd()
		for folder in folders:
			os.chdir(folder)
			os.system('sbatch sbatchMulti.bash')
		os.chdir(cwd)

	
