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
	job_set_name = "large_aggregate_test"
	folder_name_scheme = "N_"

	job_set_folder = curr_folder + 'jobs/' + job_set_name + '/'
	if not os.path.exists(job_set_folder):
		os.makedirs(job_set_folder)
	else:
		print("Job set '{}' already exists.".format(job_set_folder))

	#Make an array of what you want to vary
	attempts_for_N = [\
					[],\
					]
	Temps = [1,3,10,30,100,300,1000]
	N = [30,100,300]

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
		input_json['temp'] = Temps[index]
		# input_json['Nballs'] = N[index]
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
	# #loop folders again and run stuff
	# for job in folders:
	# 	os.system("cd {}".format(job))
	# 	os.system("pyth")

	# os.system("cd {}".format(curr_folder))

	
