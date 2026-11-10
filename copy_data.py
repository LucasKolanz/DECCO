import numpy as np
# import matplotlib.pyplot as plt
import os 
import sys
import json
import glob

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u




def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	source_prefolder = path + 'jobsCosine/lognormrelax_'
	dest_prefolder = path + 'jobs/BAPA_'
	
	
	attempts = [i for i in range(30)]

	source_dir_pattern = source_prefolder+"$a$/N_300/T_1000/"
	dest_dir_pattern = dest_prefolder+"$a$/M_1/N_300/T_1000/"


	for a,attempt in enumerate(attempts):
		source_dir = source_dir_pattern.replace("$a$",f"{attempt}")
		dest_dir = dest_dir_pattern.replace("$a$",f"{attempt}")
		if not os.path.isdir(dest_dir): 
			os.makedirs(dest_dir)
		os.system(f"cp {source_dir}* {dest_dir}.")
		print(source_dir)
		print(dest_dir) 



if __name__ == '__main__':
	main()