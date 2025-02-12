"""
This file was originally written for SpaceLab/DECCO to check the seeds of all jobs in folders matching a specified pattern.

Author: Lucas Kolanz

This file is meant to go through all jobs contained in a specified file pattern and print out their seeds. All of these seeds should
be different, and this file is meant to help verify that.

TODO: use the seeds from input.json, not seedFile.txt. Save all seeds in a list to verify they are all unique.

"""





import os
import json
import glob
import numpy as np
import subprocess
import sys
import utilities.utils as u

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u



def main():

	#Open SpaceLab default file for directory information
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	# folder_pattern_to_change = input_json["data_directory"]+"jobsCosine/lognorm_relax*"
	# new_pattern = input_json["data_directory"]+"jobsCosine/lognormrelax_{a}"


	folder_pattern_to_change = input_json["data_directory"]+"jobsCosine/lognorm*"
	new_pattern = input_json["data_directory"]+"jobsCosine/lognorm_{a}"

	


	source_folders = glob.glob(folder_pattern_to_change)

	source_folders = [i for i in source_folders if ("lognormrelax" not in i.split("/")[-1])]



	

	for folder in source_folders:
		#extract index
		index = "".join([s for s in folder.split('/')[-1] if s.isdigit()])
		
		new_folder = new_pattern.replace('{a}',index)

		# print(folder)
		# print(new_pattern.replace('{a}',index))

		os.rename(folder,new_folder)


if __name__ == '__main__':
	main()