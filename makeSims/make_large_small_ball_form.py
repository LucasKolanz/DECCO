import os
import json
import multiprocessing as mp
import subprocess

def run_job(location,num_balls): #what happens in here can be anything, not just running a subprocess
	cmd = ["python3", "{}run_sim.py".format(location), location, str(num_balls)]
	# print(cmd)
	# exit(0)
	subprocess.run(cmd)

if __name__ == '__main__': #Important for multiprocessing to put your main function in this way
	#make new output folders
	curr_folder = os.getcwd() + '/'
	job_set_name = "cuttoff_test"
	folder_name_scheme = "c_"

	try:
		subprocess.run(["make","-C","ColliderSingleCore"], check=True)
	except:
		print('compilation failed')
		exit(-1)

	job_set_folder = curr_folder + 'jobs/' + job_set_name + '/'
	if not os.path.exists(job_set_folder):
		os.makedirs(job_set_folder)
	else:
		print("Job set '{}' already exists.".format(job_set_folder))


	#Make an array of what you want to vary
	Temps = [1]
	N = [2]
	cuttoff = [1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6]
	cuttoff = [1.21,1.22,1.23,1.24,1.26,1.27,1.28,1.29]
	cuttoff = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,1.7,1.8,1.9,2.0]
	cuttoff = [5.0]
	cuttoff = [1.2]

	#actually make the folders
	folder_values = cuttoff
	folders = [job_set_folder + folder_name_scheme + str(val) + '/' for val in folder_values]
	for new_folder in folders:
		if not os.path.exists(new_folder):
			os.mkdir(new_folder)
		else:
			print("Job '{}' already exists.".format(new_folder))

	#load default input file
	with open(curr_folder+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	#loop through folders and generate input file and run_script.py for each
	#changing the necessary values in the input file
	for index,value in enumerate(folder_values):
		####################################
		######Change input values here######
		input_json['temp'] = Temps[0]
		input_json['seed'] = 100
		input_json['genBalls'] = 2
		input_json['radiiDistribution'] = 'constant'
		input_json['simTimeSeconds'] = 1e-3
		input_json['radiiFraction'] = value
		####################################
		with open(folders[index] + "input.json",'w') as fp:
			json.dump(input_json,fp,indent=4)

		#add run script and executable to folders
		os.system("cp default_files/run_sim.py {}run_sim.py".format(folders[index]))
		os.system("cp ColliderSingleCore/ColliderSingleCore.o {}ColliderSingleCore.o".format(folders[index]))

	if len(N) != len(folders):
		N = [str(N[0]) for i in range(len(folders))]

	inputs = list(zip(folders,N))
	print(inputs)

	with mp.Pool(processes=len(folders)) as pool:
		pool.starmap(run_job,list(inputs)) 


	
