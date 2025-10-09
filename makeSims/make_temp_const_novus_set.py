import os
import json
import multiprocessing as mp
import subprocess
import sys
import argparse
# import random
# import re


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

	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"Collider"], check=True)
	except:
		print('compilation failed')
		exit(-1)


	job_set_name = "constrollingfric"
	# folder_name_scheme = "T_"

	# runs_at_once = 7
	runs_at_once = 5

	attempts = [i for i in range(30)]
	# attempts = [19]



	node = 1
	N = [300]
	Temps = [3,10,30,100,300,1000]
	# Temps = [100]


	folders = []
	for n in N:
		threads = 1
		# if n == 30:
		# 	threads = 1
		# elif n == 100:
		# 	threads = 2
		# else:# n == 300:
		# 	threads = 16

		for attempt in attempts:
			for Temp in Temps:
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				job = input_json["data_directory"] + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				# job = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/fixme/"

				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				job_name = f"a={attempt},n={n},t={Temp}"


				if os.path.exists(job+"timing.txt"):
					print("Sim already complete")
				elif u.on_queue(job_name):
					print(f"Sim already on queue: {job}")
				else:
					#load default input file

					####################################
					######Change input values here######
					input_json['temp'] = Temp
					input_json['seed'] = u.rand_int()
					input_json['radiiDistribution'] = 'constant'
					input_json['N'] = n
					input_json['h_min'] = 0.5
					input_json['dataFormat'] = "h5"
					input_json['output_folder'] = job
					input_json['OMPthreads'] = threads
					input_json['impactParameter'] = -1.0
					input_json['u_s'] = 0.1
					input_json['u_r'] = 0.001
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
					sbatchfile += f"#SBATCH -J {job_name}\n"
					sbatchfile += f"#SBATCH -N {node}\n"
					sbatchfile += f"#SBATCH -n {node}\n"
					sbatchfile += f"#SBATCH -c {threads}\n\n"
					# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)


					# # sbatchfile += "#SBATCH -G {}\n".format(node)
					# # sbatchfile += 'module load gpu\n'

					# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threads)
					# sbatchfile += 'export SLURM_CPU_BIND="cores"\n'
					# # sbatchfile += 'module load hdf5/1.14.3\n'
					# sbatchfile += 'module load hdf5/1.10.8\n'
					
					# # sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
					# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"


					
					# with open(job+"sbatchMulti.bash",'w') as sfp:
					# 	sfp.write(sbatchfile)

					#add run script and executable to folders
					# os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
					os.system(f"cp {project_path}Collider/Collider.x {job}Collider.x")
					os.system(f"cp {project_path}Collider/Collider.cpp {job}Collider.cpp")
					os.system(f"cp {project_path}Collider/ball_group.hpp {job}ball_group.hpp")

					folders.append(job)


	print(folders)

	if args.run:
		cwd = os.getcwd()
		for folder in folders:
			os.chdir(folder)
			os.system('sbatch sbatchMulti.bash')
		os.chdir(cwd)








	
