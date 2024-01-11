import os
import json
import multiprocessing as mp
import subprocess

def run_job(location):
	cmd = ["python3", "{}run_sim.py".format(location), location, '7']
	# exit(0)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'
	job_set_name = "test"
	folder_name_scheme = "T_"

	attempts = [18] 
	attempts_at_once = 1
	for i in range(0,len(attempts),attempts_at_once):
		folders = []
		for j in range(i,i+attempts_at_once):
			# print(i,j)
			job_set_folder = curr_folder + 'jobs/' + job_set_name + str(attempts[j]) + '/'
			if not os.path.exists(job_set_folder):
				os.makedirs(job_set_folder)
			else:
				print("Job set '{}' already exists.".format(job_set_folder))

			#Make an array of what you want to vary
			Temps = [300]

			#actually make the folders
			folder_values = Temps
			new_folders = [job_set_folder + folder_name_scheme + str(val) + '/' for val in folder_values]
			for new_folder in new_folders:
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
				input_json['temp'] = value
				with open(new_folders[index] + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				#add run script and executable to folders
				os.system("cp default_files/run_sim.py {}run_sim.py".format(new_folders[index]))
				os.system("cp ColliderSingleCore/ColliderSingleCore.o {}ColliderSingleCore.o".format(new_folders[index]))
			folders.extend(new_folders)
		# print(folders)
		# print(len(folders))
		# exit(0)
		with mp.Pool(processes=len(folders)) as pool:
			pool.map(run_job,folders) 


	
