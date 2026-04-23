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



def delete_saved_data(directory,indices):
	if isinstance(indices, int):
		indices = [indices]
	for index in indices:
		files = glob.glob(f"{directory}{index}_*")
		[os.remove(file) for file in files]

def verify_groups_integrity(groups, M):
    _, counts = np.unique(groups, return_counts=True)
    return np.all(counts <= M)

#makes sure there aren't more then M balls per group
def verify_groups(directory,delete=False):
	indices = u.get_all_indices(directory)
	M = u.M_from_directory(directory)
	#is there some data that is possibly corrupted but doesnt have a checkpoint?
	#if so, add this to the list of indices to check
	if len(glob.glob(f"{directory}{max(indices)+M}_*")) > 0:
		indices = indices + [max(indices)+M]

	print(f"Verifying aggregation in {directory}.")
	for i in indices:
		groups = u.get_groups(directory,data_index=i,relax=False)
		if not verify_groups_integrity(groups,M):
			print(f"Groups corrupt at index {i}.")
			if delete:
				delete_saved_data(directory,i)

def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	
	data_folders = []
	data_folders = [path + 'jobs/BAPAWELD_*']


	DELETE = False


	possible_dirs = []
	for data_folder in data_folders:
		possible_dirs.extend(u.get_directores_containing(data_folder,[]))

	# possible_dirs = ["/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/BAPAWELD_9/M_3/N_300/T_1000/"]
	for d_i,directory in enumerate(possible_dirs):
		verify_groups(directory,delete = DELETE)
		


if __name__ == '__main__':
	main()
