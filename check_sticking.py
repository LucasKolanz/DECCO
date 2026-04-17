'''
This script goes through all jobs in a given folder and ensures they are connected at the end.
There is a check in the sim_looper function, but it is possible to turn it off. If that 
happens on accident, this script is handy.
'''

import sys
import os
import json
import numpy as np

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

#returns true if every ball has at least one contact
def verify_touching(pos,radii):
	# Pairwise displacement vectors: shape (n, n, d)
    diff = pos[:, None, :] - pos[None, :, :]

    # Squared pairwise distances: shape (n, n)
    dist2 = np.sum(diff * diff, axis=-1)

    # Pairwise sum of radii: shape (n, n)
    rsum = radii[:, None] + radii[None, :]

    # Touching matrix
    touching = dist2 < (rsum * rsum)

    # Ignore self-comparisons
    np.fill_diagonal(touching, False)

    # Every particle must touch at least one other
    return np.all(np.any(touching, axis=1))
	# length = radii.size
	# touching = np.zeros((length,length))

	# for i in range(0,length):
	# 	for j in range(i+1,length):
	# 		if (radii[i]+radii[j]) > np.linalg.norm(pos[i] - pos[j]): # touching
	# 			touching[i,j] = 1
	# 			touching[j,i] = 1

	# return np.sum(np.where(np.sum(touching,axis=1) > 0,1,0)) == length

def verify_aggregation(directory):
	indices = u.get_all_indices(directory)
	print(f"Verifying aggregation in {directory}")
	for i in indices:
		pos,radii,mass,moi = u.get_data(directory,data_index=i,relax=False)
		if not verify_touching(pos,radii):
			print(f"Directory not touching at index {i}\n\t{directory}")

def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	
	data_folders = []
	data_folders = [path + 'jobs/BAPA_*']


	possible_dirs = []
	for data_folder in data_folders:
		possible_dirs.extend(u.get_directores_containing(data_folder,[]))

	for d_i,directory in enumerate(possible_dirs):
		verify_aggregation(directory)
		


if __name__ == '__main__':
	main()
