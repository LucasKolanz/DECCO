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


	job_set_name = "hminTests"

	runs_at_once = 10
	attempts = [2] 
	# attempts = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	# attempts = [i for i in range(11,17)]

	# print(attempts)
	# exit(0)

	N = [2]
	h_mins = [.2,.3,.4,.5,.6,.7,.8,.9]
	folders = []
	folders_N = []
	for h_min in h_mins:
		for attempt in attempts:
			for n in N:
				job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'hmin_' + str(h_min) + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				#load default input file
				with open(curr_folder+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				####################################
				######Change input values here######
				input_json['temp'] = 1
				input_json['seed'] = 100
				input_json['genBalls'] = 2
				input_json['radiiDistribution'] = 'constant'
				# input_json['scaleBalls'] = 5e4
				# input_json['timeResolution'] = 1e2
				input_json['simTimeSeconds'] = 1e-3
				input_json['h_min'] = h_min
				# input_json['simTimeSeconds'] = 0.5e-3 # Original one
				input_json['radiiFraction'] = 2
				input_json['note'] = 'h_min of {}'.format(h_min)
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				os.system("cp default_files/run_sim.py {}run_sim.py".format(job))
				os.system("cp ColliderSingleCore/ColliderSingleCore.o {}ColliderSingleCore.o".format(job))
				folders.append(job)
				folders_N.append(n)
	# print(folders)
	# if len(N) != len(folders):
	# 	for i in range(len(folders))
	# 	N = [str(N[0]) for i in range(len(folders))]

	inputs = list(zip(folders,folders_N))
	print(inputs)

	for i in range(0,len(folders),runs_at_once):
		with mp.Pool(processes=runs_at_once) as pool:
			pool.starmap(run_job,inputs[i:i+runs_at_once]) 


	
