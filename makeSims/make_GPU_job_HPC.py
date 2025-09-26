import os
import json
import multiprocessing as mp
import subprocess
import random
import re
import sys

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
		

	# job_set_name = "TESTBAPA"
	job_set_name = "LARGEBAPAGPU1"

	# folder_name_scheme = "T_"

	SPECIAL_FOLDER = ""#"/home/lucas/Desktop/SpaceLab_data/largejob/"

	# runs_at_once = 10
	# attempts = [2] 
	attempts = [0]
	# attempts = [0]#[0,1,2,3,4,5,6,7,8,9]#,11,12,13,14,15,16,17,18,19,20] 

	N = [300*15] #final size
	M = [300] #starting sizes
	# M = [15] 
	threads = []
	# Temps = [3,10,30,100,300,1000]
	Temps = [3]
	folders = []

	totalNodes = 1
	MPITasksPerNode = 1
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 1

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_template = input_json["data_directory"] + 'jobs/' + job_set_name + '_{a}/M_{m}/N_{n}/T_{t}/'

	for attempt in attempts:
		for m in M:
			for n in N:
				for Temp in Temps:
					#load default input file
					# with open(project_path+"default_files/default_input.json",'r') as fp:
					# 	input_json = json.load(fp)
					
					# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
					job = job_template.replace('{a}',str(attempt)).replace('{m}',str(m)).replace('{n}',str(n)).replace('{t}',str(Temp))

					# print(rand_int())
					if not os.path.exists(job):
						os.makedirs(job)
					else:
						print(f"Job already exists: {job}")


					if os.path.exists(job+"timing.txt"):
						print(f"Sim already complete: {job}")

					elif u.on_queue(job):
						print(f"Sim already on queue: {job}")
					else:
						print(f"(Re)Starting job: {job}")
						####################################
						######Change input values here######
						# input_json['temp'] = Temp
						input_json['N'] = n
						input_json['M'] = m
						input_json['output_folder'] = job
						input_json['OMPthreads'] = 1#threadsPerTask
						input_json['MPInodes'] = 1
						input_json['impactParameter'] = -1.0

						# input_json['seed'] = u.rand_int()
						input_json['seed'] = 101

						# input_json['radiiDistribution'] = 'logNormal'
						# input_json['h_min'] = 0.5
						
						# input_json['timeResolution'] = 1e-6
						input_json['v_custom'] = '5'

						# input_json['simTimeSeconds'] = 1e-6
						input_json['simTimeSeconds'] = 5e-4

						input_json['dataFormat'] = "csv"
						input_json['simType'] = "BAPA"
						input_json['random_folder_template'] = input_json['data_directory']+"jobs/const300/"

						# input_json['u_s'] = 0.5
						# input_json['u_r'] = 0.5
						# input_json['note'] = "Does this work at all?"
						####################################

						with open(job + "input.json",'w') as fp:
							json.dump(input_json,fp,indent=4)


						sbatchfile = ""
						sbatchfile += "#!/bin/bash\n"
						# sbatchfile += "#SBATCH -C gpu\n"
						# sbatchfile += "#SBATCH -q regular\n"
						# sbatchfile += "#SBATCH -t 0:10:00\n"
						# sbatchfile += f'#SBATCH --partition=dri.q\n'

						# #FOR ENGR CLUSTER
						# #NAME ORDER needs to be same as the file path order
						# sbatchfile += f"#SBATCH -J LBGPU,GPUS={totalMPITasks},a={attempt},m={m},n={n},t={Temp}\n"
						# # sbatchfile += "#SBATCH -A kolanzl\n"
						# # sbatchfile += "#SBATCH --partition=share\n"
						# sbatchfile += "#SBATCH --partition=dgxh\n"
						# sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
						# sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
						# sbatchfile += f"#SBATCH --gpus-per-node=1\n"
						# # sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
						# # sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

						# # sbatchfile += "#SBATCH -G {}\n".format(node)
						# # sbatchfile += 'module load gpu\n'

						# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
						# sbatchfile += 'lscpu\n'
						# # sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
						# # sbatchfile += 'module load hdf5/1.14.3\n'
						# # sbatchfile += 'module load hdf5/1.10.8\n'
						# # sbatchfile += 'module load gcc/8.3\n'
						# # sbatchfile += 'module load gnu12/12.3.0\n'

						# sbatchfile += 'module purge\n'
						# sbatchfile += 'module load hdf5/1.10.5_mpich-3.3\n'
						# sbatchfile += 'module load slurm/24.05\n'
						# sbatchfile += 'module load python/3.11\n'
						# sbatchfile += 'module load nvhpcsdk/2024\n'
						# sbatchfile += 'module load nvhpc/25.3\n'
						# sbatchfile += 'module load mpich/3.3\n'
						# sbatchfile += 'module list\n'
						# sbatchfile += 'echo "Running on host: $(hostname)"\n'
						# sbatchfile += 'export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH\n'
						# sbatchfile += 'echo $LD_LIBRARY_PATH\n'
						# sbatchfile += 'ls /usr/lib64/libc.*\n'
						# sbatchfile += "unset NV_ACC_DEBUG\n"
						# sbatchfile += "unset NV_ACC_NOTIFY\n"
						# # sbatchfile += "export NV_ACC_NOTIFY=3\n"
						# # sbatchfile += "export NV_ACC_DEBUG=0x800\n"
						# sbatchfile += "export ACC_DEVICE_TYPE=nvidia\n"
						# # sbatchfile += "export CUDA_LAUNCH_BLOCKING=1\n"
						# sbatchfile += "nvidia-smi\n"

						# #NAME ORDER needs to be same as the file path order
						# sbatchfile += f"#SBATCH -J LBGPU,GPUS={totalMPITasks},a={attempt},m={m},n={n},t={Temp}\n"
						# # sbatchfile += "#SBATCH -A kolanzl\n"
						# # sbatchfile += "#SBATCH --partition=share\n"
						# sbatchfile += "#SBATCH --partition=dgxh\n"
						# sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
						# sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
						# sbatchfile += f"#SBATCH --gpus-per-node=1\n"
						# # sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
						# # sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

						# # sbatchfile += "#SBATCH -G {}\n".format(node)
						# # sbatchfile += 'module load gpu\n'

						# sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
						# sbatchfile += 'lscpu\n'
						# os.system(f"mkdir {job}/lib")
						# os.system(f"cp /usr/lib64/libhdf5.so.103 {job}/lib/.")
						# os.system(f"cp /usr/lib64/libhdf5_cpp.so.103 {job}/lib/.")


						#NOVUS CLUSTER
						sbatchfile = ""
						sbatchfile += "#!/bin/bash\n"
						# sbatchfile += "#SBATCH -C gpu\n"
						# sbatchfile += "#SBATCH -q regular\n"
						# sbatchfile += "#SBATCH -t 0:10:00\n"
						# sbatchfile += f'#SBATCH --partition=dri.q\n'

						#FOR ENGR CLUSTER
						#NAME ORDER needs to be same as the file path order
						sbatchfile += f"#SBATCH -J LBGPU,GPUS={totalMPITasks},a={attempt},m={m},n={n},t={Temp}\n"
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

						# os.system(f"cp /usr/lib64/libdl.so.2 {job}/lib/.")
						# os.system(f"cp /usr/lib64/libsz.so.2 {job}/lib/.")
						# os.system(f"cp /usr/lib64/libz.so.1 {job}/lib/.")
						# os.system(f"cp /usr/lib64/libzma.so.5 {job}/lib/.")



						randint = random.randint(0, 29)
						# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm_relax{randint}/N_30/T_3/27_RELAXconstants.csv {job}{m}_constants.csv")
						# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm_relax{randint}/N_30/T_3/27_RELAXsimData.csv {job}{m}_simData.csv")
						# os.system(f"cp /media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm_relax{randint}/N_30/T_3/27_RELAXenergy.csv {job}{m}_energy.csv")
						# source = "/media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm{randint}/N_30/T_3/{m}_*"
						# if M == 3:
							# source = "/media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm{randint}/N_30/T_3/2_R*"
						if not os.path.exists(f"{job}{m}_RELAXdata.h5"):
							os.system(f"cp {input_json['data_directory']}jobs/const300/{m}_data.h5 {job}{m}_data.h5")
							# os.system(f"cp {input_json['data_directory']}jobs/const300/{m}_constants.csv {job}{m}_constants.csv")
							# os.system(f"cp {input_json['data_directory']}jobs/const300/{m}_simData.csv {job}{m}_simData.csv")
							# os.system(f"cp {input_json['data_directory']}jobs/const300/{m}_energy.csv {job}{m}_energy.csv")
							
						folders.append(job)
	# print(folders)


	print(folders)
	cwd = os.getcwd()
	for folder in folders:
		os.chdir(folder)
		os.system('sbatch sbatchMulti.bash')
	os.chdir(cwd)

	
