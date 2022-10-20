import os
import json
import multiprocessing as mp
import subprocess

def run_job(location,num_balls):
	cmd = ["python3", "{}run_sim_parallel.py".format(location), location, str(num_balls)]
	# print(cmd)
	# exit(0)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	#switch directories and make the c++ file
	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","ColliderParallel"], check=True)
	except:
		print('compilation failed')
		exit(-1)

	job_set_name = "parallel_test"
	folder_name_scheme = "N_"

	job_set_folder = curr_folder + 'jobs/' + job_set_name + '/'
	if not os.path.exists(job_set_folder):
		os.makedirs(job_set_folder)
	else:
		print("Job set '{}' already exists.".format(job_set_folder))

	#Make an array of what you want to vary
	Temps = [1000]
	N = [1000]

	#actually make the folders
	folder_values = N
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
		###########enter input values here###########
		input_json['temp'] = Temps[index]
		# input_json['gridSize'] = 1.0
		# input_json['Nballs'] = N[index]
		#############################################
		with open(folders[index] + "input.json",'w') as fp:
			json.dump(input_json,fp,indent=4)

		#add run script and executable to folders
		os.system("cp default_files/run_sim_parallel.py {}run_sim_parallel.py".format(folders[index]))
		os.system("cp ColliderParallel/ColliderSingletoParallel.o {}ColliderSingletoParallel.o".format(folders[index]))

	if len(N) != len(folders):
		N = [str(N[0]) for i in range(len(folders))]

	inputs = list(zip(folders,N))
	print(inputs)

	with mp.Pool(processes=len(folders)) as pool:
		pool.starmap(run_job,list(inputs)) 
	# #loop folders again and run stuff
	# for job in folders:
	# 	os.system("cd {}".format(job))
	# 	os.system("pyth")

	# os.system("cd {}".format(curr_folder))

	
