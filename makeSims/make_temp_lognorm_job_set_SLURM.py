import os
import json
import multiprocessing as mp
import subprocess

relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def run_job(location,num_balls):
	cmd = ["python3", "{}run_sim.py".format(location), location, str(num_balls)]
	# print(cmd)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C",project_path+"ColliderSingleCore"], check=True)
	except:
		print('compilation failed')
		exit(-1)


	job_set_name = "lognorm"
	# folder_name_scheme = "T_"

	# runs_at_once = 7
	# attempts = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [i for i in range(10)]
	# attempts = [1] 
	attempts_300 = [i for i in range(5)]

	#test it out first
	# attempts = [0]
	# attempts_300 = [0]

	N = [300,30,100]
	# N = [5]
	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	folders = []
	for n in N:
		for Temp in Temps:
			temp_attempt = attempts
			if n == 300:
				temp_attempt = attempts_300
			for attempt in temp_attempt:
				job = project_path + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				#load default input file
				with open(project_path+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['seed'] = 'default'
				input_json['radiiDistribution'] = 'logNormal'
				input_json['h_min'] = 0.5
				# input_json['u_s'] = 0.5
				# input_json['u_r'] = 0.5
				input_json['note'] = "Runs testing h_min = 0.5 (5e-6) with lognormal distribution"
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)


				sbatchfile = ""
				sbatchfile += "#!/bin/bash\n"
				sbatchfile += "#SBATCH -A m2651\n"
				sbatchfile += "#SBATCH -C gpu\n"
				sbatchfile += "#SBATCH -q regular\n"
				sbatchfile += "#SBATCH -t 0:10:00\n"
				sbatchfile += "#SBATCH -J {}\n".format(job_set_name)
				sbatchfile += "#SBATCH -N {}\n".format(1)#(node)
				# sbatchfile += "#SBATCH -G {}\n".format(node)
				# sbatchfile += "#SBATCH -c {}\n\n".foramt(2*thread)
				# sbatchfile += 'module load gpu\n'
				# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(thread)
				sbatchfile += 'export SLURM_CPU_BIND="cores"\n'
				
				# sbatchfile += "srun -n {} -c {} --cpu-bind=cores numactl --interleave=all ./ColliderMultiCore.x {} 2>sim_err.log 1>sim_out.log".format(node,thread*2,job)
				sbatchfile += "srun -n {} -c {} ./ColliderSingleCore.o {} {} 2>sim_err.log 1>sim_out.log".format(node,thread*2,job,n)


				
				with open(job+"sbatchMulti.bash",'w') as sfp:
					sfp.write(sbatchfile)

				#add run script and executable to folders
				os.system(f"cp {project_path}default_files/run_sim.py {job}run_sim.py")
				os.system(f"cp {project_path}ColliderSingleCore/ColliderSingleCore.x {job}ColliderSingleCore.x")
				os.system(f"cp {project_path}ColliderSingleCore/ColliderSingleCore.cpp {job}ColliderSingleCore.cpp")
				# os.system("cp default_files/run_multicore_sim.py {}run_multicore_sim.py".format(job))
				# os.system("cp ColliderMultiCore/ColliderMultiCore.x {}ColliderMultiCore.x".format(job))
				os.system(f"cp {project_path}ColliderMultiCore/ball_group_multi_core.hpp {job}ball_group_multi_core.hpp")
				# if input_json['simType'] != "BPCA":
				# 	os.system("cp ../jobs/collidable_aggregate_1200/* {}".format(job))

				folders.append(job)

print(folders)
# cwd = os.getcwd()
# for folder in folders:
# 	os.chdir(folder)
# 	os.system('sbatch sbatchMulti.bash')
# os.chdir(cwd)







	
