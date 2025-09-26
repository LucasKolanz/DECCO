import os
import json
import multiprocessing as mp
import subprocess

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
		
	job_set_name = "videoCollision"

	# folder_name_scheme = "T_"

	# SPECIAL_FOLDER = ""#"/home/lucas/Desktop/SpaceLab_data/largejob/"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [6] 
	folders = []

	totalNodes = 1
	MPITasksPerNode = 2
	totalMPITasks = totalNodes*MPITasksPerNode
	threadsPerTask = 24
	# threadsPerTask = 48
	for attempt in attempts:

		#load default input file
		with open(project_path+"default_files/default_input.json",'r') as fp:
			input_json = json.load(fp)
		
		# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
		job = input_json["data_directory"] + 'jobs/' + job_set_name + '_' + str(attempt) + '/'
		
		if not os.path.exists(job):
			os.makedirs(job)
		else:
			print("Job '{}' already exists.".format(job))

		####################################
		######Change input values here######
		# input_json['temp'] = 3
		# input_json['N'] = 5
		input_json['output_folder'] = job
		input_json['OMPthreads'] = threadsPerTask
		input_json['MPInodes'] = totalMPITasks

		# input_json['simTimeSeconds'] = 0.001  
		# input_json['timeResolution'] = 1e-5

		input_json['simTimeSeconds'] = 0.001
		input_json['timeResolution'] = (1e-5)/5.0

		input_json['seed'] = 100
		input_json['v_custom'] = 100
		# input_json['h_min'] = 0.5
		input_json['dataFormat'] = "csv"
		# input_json['simType'] = "BPCA"
		input_json['simType'] = "collider"

		# input_json['projectileName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP/297_simData.csv"
		# input_json['targetName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP/297_simData.csv"
		input_json['projectileName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP/297_RELAXsimData.csv"
		input_json['targetName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP/297_RELAXsimData.csv"
		# # input_json['projectileName'] = "/media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm_relax7/N_300/T_3/297_RELAXsimData.csv"
		# input_json['targetName'] = "/media/kolanzl/easystore/SpaceLab_data/jobsCosine/lognorm_relax7/N_300/T_3/297_RELAXsimData.csv"

		# input_json['projectileName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP1/27_RELAXsimData.csv"
		# input_json['targetName'] = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/colliderTP1/27_RELAXsimData.csv"
		

		# input_json['u_s'] = 0.5
		# input_json['u_r'] = 0.5
		# input_json['note'] = "Does this work at all?"
		####################################

		with open(job + "input.json",'w') as fp:
			json.dump(input_json,fp,indent=4)

		sbatchfile = ""
		sbatchfile += "#!/bin/bash\n"
		# sbatchfile += "#SBATCH -A m2651\n"
		# sbatchfile += "#SBATCH -C gpu\n"
		# sbatchfile += "#SBATCH -q regular\n"
		# sbatchfile += "#SBATCH -t 0:10:00\n"
		# sbatchfile += f'#SBATCH --partition=dri.q\n'
		sbatchfile += f"#SBATCH -J coll\n"
		sbatchfile += f"#SBATCH --nodes {totalNodes}\n"
		sbatchfile += f"#SBATCH --ntasks-per-node {totalMPITasks}\n"
		sbatchfile += f"#SBATCH --cpus-per-task {threadsPerTask}\n\n"
		# sbatchfile += "#SBATCH -N {}\n".format(1)#(node)

		# sbatchfile += "#SBATCH -G {}\n".format(node)
		# sbatchfile += 'module load gpu\n'

		sbatchfile += 'export OMP_NUM_THREADS={}\n'.format(threadsPerTask)
		# sbatchfile += 'export SLURM_CPU_BIND="socket"\n'
		# sbatchfile += 'module load hdf5/1.14.3\n'
		# sbatchfile += 'module load hdf5/1.10.8\n'
		sbatchfile += 'module load gnu12/12.3.0\n'
		sbatchfile += 'module load openmpi4/4.1.6\n'
		# sbatchfile += 'module swap openmpi4/4.1.6 mpich\n'

		
		# sbatchfile += f"srun -n {node} -c {threads} --cpu-bind=cores numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
		# sbatchfile += f"srun --ntasks-per-node={MPITasksPerNode} --cpus-per-task={threadsPerTask} --cpu-bind=socket numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
		# sbatchfile += f"mpirun --bind-to socket --map-by node numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
		# sbatchfile += f"mpirun -n {totalMPITasks} numactl --interleave=all {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"
		sbatchfile += f"mpirun -n {totalMPITasks} {job}Collider.x {job} 2>>sim_err.log 1>>sim_out.log\n"


		
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

		# os.system(f"cp /home/lucas/Desktop/SpaceLab_data/test2/N_5/T_3/*data.h5 {job}.")
		
		
	# print(folders)





###############For running locally
	# print(folders)

	# # for i in range(0,len(folders),runs_at_once):
	# # 	with mp.Pool(processes=runs_at_once) as pool:
	# # 		pool.starmap(run_job,inputs[i:i+runs_at_once]) 
	
	# with mp.Pool(processes=runs_at_once) as pool:
	# 	for folder in folders:
	# 		# input_data = inputs[i:i+runs_at_once]
	# 		pool.apply_async(run_job, (folder,))

	# 	pool.close()
	# 	pool.join()

	# print(folders)
	# cwd = os.getcwd()
	# for folder in folders:
	# 	os.chdir(folder)
	# 	os.system('qsub qsub.bash')
	# os.chdir(cwd)

	
