import os
import json
import multiprocessing as mp
import subprocess

def run_job(location,num_balls):
	cmd = ["python3", "{}run_sim_parallel.py".format(location), location, str(num_balls)]
	# print(cmd)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","ColliderParallel"], check=True)
	except:
		print('compilation failed')
		exit(-1)

	# job_set_name = "LargePairParallelTest"
	job_set_name = "lognorm_radius_test"
	job_set_name = "test"
	
	# folder_name_scheme = "T_"
	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [1] 
	N = [10]
	Temps = [1]
	folders = []
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				#load default input file
				with open(curr_folder+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				####################################
				######Change input values here######
				input_json['temp'] = Temp
				# input_json['gridSize'] = 4e-5
				input_json['seed'] = 100
				input_json['radiiDistribution'] = 'logNormal'
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				os.system("cp default_files/run_sim_parallel.py {}run_sim_parallel.py".format(job))
				os.system("cp ColliderParallel/ColliderParallel.o {}ColliderParallel.o".format(job))
				folders.append(job)
	# print(folders)
	if len(N) != len(folders):
		N = [str(N[0]) for i in range(len(folders))]

	inputs = list(zip(folders,N))
	print(inputs)

	for i in range(0,len(folders),runs_at_once):
		with mp.Pool(processes=runs_at_once) as pool:
			pool.starmap(run_job,inputs[i:i+runs_at_once]) 


	
