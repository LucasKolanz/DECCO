##restart jobs that have errored

#Error 1: restart failed, number of balls didn't carry over so there are an incorrect 
#			number of colums in at least 1 sim_data output
#Error 2: sim fails to write all the data (so there aren't the correct number of colums in simData)
#Error *: possible signed integer over flow in number of steps for a sim
#			Indicated by a specific sequence at the end of sim_errors.txt
#Error general: did we get "Simulation complete!" within the last 10 lines of sim_error.log

import os
import glob
import numpy as np
import subprocess
import check_for_errors as cfe
import check_for_errors_h5 as cfeh5

# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def modify_last_line(file_path, new_content):
    with open(file_path, 'r+') as file:
        lines = file.readlines()
        if lines:
            lines[-1] = new_content
        
        file.seek(0)
        file.writelines(lines)
        file.truncate()

def restart_job(folder,test=True,move_folder=''):
	if test:
		print('==================================================================================')
		print("IN TEST MODE")
		print('==================================================================================')
	if len(move_folder) > 0: #move data to new folder specified in move_folder
		if os.path.exists(move_folder): #if move_folder already exists
			#if it already exists then we need to change the name of it so it doesn't overwrite
			move_folder = move_folder[:-1] + '_MOVE-'
			move_index = 0
			while os.path.exists(move_folder+str(move_index)): #check for lowest number transpher. 
				move_index += 1

			move_folder += str(move_index) + '/'
			if test:
				print(f"make {move_folder}")
			else:
				os.makedirs(move_folder)
		
		else : #make it
			if test:
				print(f"make {move_folder}")
			else:
				os.makedirs(move_folder)

		command = f"mv {folder}* {move_folder}."
		if test:
			print(command)
		else:
			os.system(command)

	else: #remove and restart
		if test:
			command = "ls"
		else:
			command = "rm"

		os.system(f"{command} {folder}*.csv")
		os.system(f"{command} {folder}*.txt")
		os.system(f"{command} {folder}*.x")
		os.system(f"{command} {folder}*.py")
		os.system(f"{command} {folder}*.cpp")
		os.system(f"{command} {folder}*.hpp")
		os.system(f"{command} {folder}*.h5")

		try:
			# os.chdir("{}ColliderSingleCore".format(curr_folder))
			subprocess.run(["make","-C","Collider"], check=True)
		except Exception as e:
			print('compilation failed')
			print(e)
			exit(-1)

		if not test:
			# os.system("cp default_files/run_sim.py {}run_sim.py".format(folder))
			os.system("cp Collider/Collider.x {}Collider.x".format(folder))
			os.system("cp Collider/Collider.cpp {}Collider.cpp".format(folder))
			os.system("cp Collider/ball_group.hpp {}ball_group.hpp".format(folder))


		cwd = os.getcwd()
		os.chdir(folder)
		if not test:
			modify_last_line("qsub.bash",f"{folder}Collider.x {folder}\n")
			os.system('qsub qsub.bash')
		else:
			os.system('ls qsub.bash')
			os.system('ls input.json')
		os.chdir(cwd)

def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	curr_folder = os.getcwd() + '/'

	job_folder = 'jobs/'###FOR RUNNIN ON COSINE/NOVUS
	job_folder = 'jobsCosine/'##FOR RUNNING ON LOCAL TESTING COSINE JOBS
	job_folder = 'jobsNovus/'##FOR RUNNING ON LOCAL TESTING NOVUS JOBS
	move_job_folder = 'erroredJobs/' ##either way move here

	job = input_json["data_directory"] + job_folder + 'lognorm$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + job_folder + 'constant$a$/N_$n$/T_$t$/'
	# move_folder = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	# attempts = [1]

	N = [30,100,300]
	# N=[300]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	if job_folder == "jobsCosine/":
		error_folders = cfe.check_error(job,cfe.errorn1,N,Temps,attempts)
		error_folders.extend(cfe.check_error(job,cfe.error0,N,Temps,attempts))
		error_folders.extend(cfe.check_error(job,cfe.error1,N,Temps,attempts))
		error_folders.extend(cfe.check_error(job,cfe.error2,N,Temps,attempts))
		error_folders.extend(cfe.check_error(job,cfe.error3,N,Temps,attempts))
		error_folders.extend(cfe.check_error(job,cfe.error4,N,Temps,attempts))
		error_folders = list(set(error_folders))
		# print(error_folders)

	elif job_folder == "jobsNovus/":
		error_folders = cfeh5.check_error(job,cfeh5.errorn1,N,Temps,attempts)
		error_folders.extend(cfeh5.check_error(job,cfeh5.error0,N,Temps,attempts))
		error_folders.extend(cfeh5.check_error(job,cfeh5.error1,N,Temps,attempts))
		error_folders.extend(cfeh5.check_error(job,cfeh5.error2,N,Temps,attempts))
		error_folders.extend(cfeh5.check_error(job,cfeh5.error3,N,Temps,attempts))
		error_folders.extend(cfeh5.check_error(job,cfeh5.error4,N,Temps,attempts))
		error_folders = list(set(error_folders))
		# print(error_folders)
	else:
		print("check job_folder in restart_error_jobs.py")

	for folder in error_folders:
		# print(folder)

		restart_job(folder,test=True,move_folder=folder.replace(job_folder,move_job_folder)) #if move_folder is specified it will move the errored jobs
		# restart_job(folder,test=False,move_folder='') #keep move_folder empty if Deleting and restarting jobs

if __name__ == '__main__':
	main()