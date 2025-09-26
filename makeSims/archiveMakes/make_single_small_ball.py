import os
import json
import multiprocessing as mp
import subprocess

def run_job(location,num_balls):
	cmd = ["python3", "{}run_sim.py".format(location), location, str(num_balls)]
	# print(cmd)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","ColliderSingleCore"], check=True)
	except:
		print('compilation failed')
		exit(-1)
		
	job_set_name = "lognorm_radius_test"
	job_set_name = "test"
	job_set_name = "small_balls"
	# folder_name_scheme = "T_"

	runs_at_once = 1
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [1] 
	N = [2]
	Temps = [1]
	folders = []
	for attempt in attempts:
		for n in N:
			for Temp in Temps:
				job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'#\
							# + 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
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
				input_json['seed'] = 100
				input_json['genBalls'] = 2
				input_json['radiiDistribution'] = 'constant'
				# input_json['scaleBalls'] = 5e4
				# input_json['timeResolution'] = 1e2
				input_json['simTimeSeconds'] = 1e-3
				input_json['h_min'] = 0.1
				# input_json['simTimeSeconds'] = 0.5e-3 # Original one
				input_json['radiiFraction'] = 2
				input_json['note'] = 'Normal h_min of 0.1'
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				os.system("cp default_files/run_sim.py {}run_sim.py".format(job))
				os.system("cp ColliderSingleCore/ColliderSingleCore.x {}ColliderSingleCore.x".format(job))
				folders.append(job)
	# print(folders)
	if len(N) != len(folders):
		N = [str(N[0]) for i in range(len(folders))]

	inputs = list(zip(folders,N))
	print(inputs)

	for i in range(0,len(folders),runs_at_once):
		with mp.Pool(processes=runs_at_once) as pool:
			pool.starmap(run_job,inputs[i:i+runs_at_once]) 


	
