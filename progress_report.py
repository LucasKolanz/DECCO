##Print out progress report of what jobs have started/finished/initialized

import os
import glob
import numpy as np
import h5py
import json
from itertools import product

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

patterns = ["$a$","$m$","$n$","$t$"]

def unroll(*lists):
    normalized = [lst if len(lst) > 0 else [None] for lst in lists]
    return list(product(*normalized))


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
			if filename.endswith(".h5") or filename.endswith("simData.csv"):
				return 1

	return 0


#gets the status of every job following the pattern given in job_base
#and the unrolled variables
def get_status(job_base,unrolled):

	stati = np.full((len(unrolled)),fill_value=np.nan);

	unfinished_jobs = []

	for a_i,attrs in enumerate(unrolled):
		job = job_base
		for v_i,value in enumerate(attrs):
			if value is not None:
				job = job.replace(patterns[v_i],str(value))
		job_status = status(job)
		if job_status == -1:
			# print(f"job doesnt exist: {job}")
			unfinished_jobs.append(job)
			# pass
		elif job_status == 0:
			# print(f"job is initialized: {job}")
			unfinished_jobs.append(job)
			# pass
		elif job_status == 1:
			# print(f"job is started: {job}")
			
			unfinished_jobs.append(job)

		stati[a_i] = job_status

	return stati, unfinished_jobs

def get_progress_bar(stati,length):

	# print(stati)
	
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

	
	# job = input_json["data_directory"] + 'jobsCosine/lognorm_relax$a$/N_$n$/T_$t$/'
	# job = input_json["data_directory"] + 'jobsCosine/lognorm$a$/N_$n$/T_$t$/'
	# job = input_json["data_directory"] + 'jobsNovus/constantX_relax$a$/N_$n$/T_$t$/'
	
	# job = input_json["data_directory"] + 'jobs/BAPA_$a$/M_$m$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobs/constrollingfricrelax$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobs/constrollingfric$a$/N_$n$/T_$t$/'
	print(job)


	attempts = [i for i in range(30)]
	# attempts = [0]

	# N = [30,100,300]
	N = [300]
	# M = [3,5,10,15]
	M=[]
	Temps = [3,10,30,100,300,1000]
	# Temps = [1000]



	#attempts needs to be first, otherwise follow the order in the list of patterns (called patterns)
	unrolled = unroll(attempts,M,N,Temps)

	job_status,unfinished_jobs = get_status(job,unrolled)

	# progress_bar_percents = np.full((len(M),len(N),len(Temps),3),fill_value=-1,dtype=np.float64)
	lines = []

	max_len = -1

	skipsize = int(len(job_status)/len(attempts))

	for start in range(skipsize):
		m = unrolled[start][1]
		n = unrolled[start][2]
		T = unrolled[start][3]
		
		length = min(len(attempts),100)
		progress_bar = get_progress_bar(job_status[start::skipsize],length)

		line_start = "Progress of "
		if m is not None:
			line_start += f"M={m}, "
		if n is not None: 
			line_start += f"N={n}, "
		if T is not None: 
			line_start += f"T={T}, "
		line = f"{line_start} [{''.join(progress_bar)}]"
		lines.append(line)
		if len(line) > max_len:
			max_len = len(line)

	# for m_i,m in enumerate(M):		
		# for n_i,n in enumerate(N):
		# 	for T_i,T in enumerate(Temps):
				
		# 		length = min(len(attempts),100)
		# 		progress_bar = get_progress_bar(job_status[:,m_i,n_i,T_i],length)
		# 		line = f"Progress of M={m}, N={n}, T={T}: [{''.join(progress_bar)}]"
		# 		lines.append(line)
		# 		if len(line) > max_len:
		# 			max_len = len(line)

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

	print(unfinished_jobs)

	


if __name__ == '__main__':
	main()
