##Print out progress report of what jobs have started/finished/initialized

import os
import glob
import numpy as np
# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'



#If job folder doesnt exist return -1, 
#	is initialized return 0, 
#	if it is started return 1, 
#	if it is finished return 2
def status(fullpath):
	if not os.path.exists(fullpath):
		return -1
	elif os.path.exists(fullpath+"timing.txt"):
		return 2
	else:
		# Loop through all files in the directory fullpath
		for filename in os.listdir(fullpath):
			if filename.endswith(".h5"):
				return 1

	return 0



def get_status(job_base,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts

	stati = np.full((len(num_attempts),len(N),len(Temps)),fill_value=np.nan);

	for a_i,attempt in enumerate(attempts):
		for n_i,n in enumerate(N):
			for T_i,Temp in enumerate(Temps):
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				stati[a_i,n_i,T_i] = status(job)

	return stati

def get_progress_bar(stati,length):
	
	where_finished = np.where(stati==2,1,0)
	where_started = np.where(stati==1,1,0)
	where_initialized = np.where(stati==0,1,0)
	where_notinitialized = np.where(stati==-1,1,0)
	percent_finished = (1.0*np.sum(where_finished))/(1.0*len(stati))
	percent_started = (1.0*np.sum(where_started))/(1.0*len(stati))
	percent_initialized = (1.0*np.sum(where_initialized))/(1.0*len(stati))
	percent_notinitialized = (1.0*np.sum(where_notinitialized))/(1.0*len(stati))

	percents = [percent_finished,percent_started,percent_initialized,percent_notinitialized]

	number_chars = []
	inc_ind = 0
	min_diff = -1
	for i,percent in enumerate(percents):
		num = np.floor(percent*length)
		diff = percent*length - np.ceil(percent*length)
		if diff < min_diff:
			min_diff = diff
			inc_ind = i
		number_chars.append(int(num))
	if np.sum(number_chars) < length:
		number_chars[inc_ind] += 1

	return number_chars[0]*"#" + number_chars[1]*"-" + number_chars[2]*"_" + number_chars[3]*" "

def main():

	curr_folder = os.getcwd() + '/'

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job = curr_folder + 'jobs/tempVarianceRand_attempt$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/weakseed$a$/N_$n$/T_$t$/'
	job = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobsNovus/testError$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsNovus/const$a$/N_$n$/T_$t$/'
	print(job)


	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[5]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	job_status = get_status(job,N,Temps,attempts)
	progress_bar_percents = np.full((len(N),len(Temps),3),fill_value=-1,dtype=np.float64)
	lines = []

	max_len = -1
	for n_i,n in enumerate(N):
		for T_i,T in enumerate(Temps):	
			
			length = min(len(attempts),100)
			progress_bar = get_progress_bar(job_status[:,n_i,T_i],length)
			line = f"Progress of N={n}, T={T}: [{''.join(progress_bar)}]"
			lines.append(line)
			if len(line) > max_len:
				max_len = len(line)

	output_line = list("_"*(max_len+2))
	output_line[0] = " "
	output_line[-1] = " "
	output = "".join(output_line) + '\n'

	for line in lines:
		output_line = list(" "*(max_len+2))
		output_line[0] = "|"
		output_line[-1] = "|"
		output += "".join(output_line) + '\n'
		
		output_line = list(" "*(max_len+2))
		output_line[0] = "|"
		output_line[-1] = "|"

		start = (max_len+2)-1-len(line)
		output_line[start:-1] = line
		output += "".join(output_line) + '\n'

		output_line = list("_"*(max_len+2))
		output_line[0] = "|"
		output_line[-1] = "|"
		output += "".join(output_line) + '\n'


	

	print(output)

	# linestart = int(len(line)/2 - len(section)/2)
	# line[linestart:len(section)] = section[:]

	# print("".join(line))

	# separator = list(totlen*"_")
	# separator[0] = "|"
	# separator[-1] = "|"
	# print("".join(separator))
	# out_lines = []
	# totlen = 50
	# top = list(totlen*"_")
	# top[0] = " "
	# top[-1] = " "
	# print("".join(top))
	


if __name__ == '__main__':
	main()