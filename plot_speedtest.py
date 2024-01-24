##Check output of simulations for errors

#Error 1: restart failed, number of balls didn't carry over so there are an incorrect 
#			number of colums in at least 1 sim_data output
#Error 2: sim fails to write all the data (so there aren't the correct number of colums in simData)
#Error *: possible signed integer over flow in number of steps for a sim
#			Indicated by a specific sequence at the end of sim_errors.txt
#Error general: did we get "Simulation complete!" within the last 10 lines of sim_error.log

import os
import glob
import numpy as np
# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


#If fullpath has error1 in it, return 1, if not return 0
def get_timing(fullpath):
	if os.path.exists(fullpath+"timing.txt"):
		with open(fullpath+"timing.txt") as fp:
			file = fp.readlines()
			time = file[-2]
		time = float(time.split(" ")[-2]) 
		return time
	else:
		return np.nan

	
	
def get_all_times(job_base,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts
	times = np.full((len(Temps),len(attempts),len(N)),fill_value=np.nan)
	# N=[30]
	# Temps = [3]
	# valid_count = 0
	for T_i,Temp in enumerate(Temps):
		for a_i,attempt in enumerate(attempts):
			for n_i,n in enumerate(N):
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				# if os.path.exists(job+"timing.txt"):
				# 	valid_count += 1
				times[T_i,a_i,n_i] = get_timing(job)
			

	return times




def main():

	curr_folder = os.getcwd() + '/'
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)


	job = input_json["data_directory"] + 'speedtests/speedtest$a$/N_$n$/T_$t$/'


	attempts = [i for i in range(30)]
	attempts = [1,2,4,8,16]

	N = [100,150,200,219]
	# N=[5]

	Temps = [3]
	# Temps = [1000]


	times = get_all_times(job,N,Temps,attempts)
	print(times)



if __name__ == '__main__':
	main()