'''
This script goes through all jobs in a given folder and ensures they are connected at the end.
There is a check in the sim_looper function, but it is possible to turn it off. If that 
happens on accident, this script is handy.
'''

import sys
import os
import json
import numpy as np
import glob

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u




#makes sure the data for N=300 exists and timing.txt
#if not, delete timing.txt and job_data.csv
def verify_finished(directory,delete=False):
	N = u.N_from_directory(directory)

	print(f"Verifying finished in {directory}.")

	if os.path.exists(directory+f"timing.txt"):
		if not os.path.exists(directory+f"{N}_checkpoint.txt"):
			print(f"\ttiming.txt exits but {N}_checkpoint.txt does not.")
			if delete:
				os.remove(directory+f"timing.txt")
				os.remove(directory+f"{N}_checkpoint.txt")
				os.remove(directory+f"job_data.csv")

def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	
	data_folders = []
	data_folders = [path + 'jobs/BAPA_*']

	DELETE = False

	possible_dirs = []
	for data_folder in data_folders:
		possible_dirs.extend(u.get_directores_containing(data_folder,[]))

	# possible_dirs = ["/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/BAPAWELD_9/M_3/N_300/T_1000/"]
	for d_i,directory in enumerate(possible_dirs):
		verify_finished(directory,delete = DELETE)
		


if __name__ == '__main__':
	main()
