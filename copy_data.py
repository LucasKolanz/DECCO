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

	source_prefolder = path + 'jobsCosine/lognorm_relax'
	dest_prefolder = path + 'jobs/BAPA_'

	# dataset_name = data_prefolder.split("/")[-1]

	# sav = path+'data/{}_averageData.csv'.format(dataset_name)
	# # figure_folder = 'figuresCompare/'
	# figure_folder = path+'data/figures/'


	# temps = [3,10,30,100,300,1000]
	# # temps = [1000]
	# Nums = [30,100,300]
	# Nums = [300]
	# M = [1]
	
	
	attempts = [i for i in range(20)]

	source_dir_pattern = source_prefolder+"$a$/N_300/T_1000/"
	dest_dir_pattern = dest_prefolder+"$a$/M_1/N_300/T_1000/"

	# for a,attempt in enumerate(attempts):

	# 	source_dir = source_dir_pattern.replace("$a$",f"{attempt}")
	# 	dest_dir = dest_dir_pattern.replace("$a$",f"{attempt}")
	# 	os.system(f"cp {source_dir}* {dest_dir}.")
	# 	# os.makedirs(dest_dir)
	# 	# print(source_dir)
	# 	# print(dest_dir) 

	# 	os


				
	for a,attempt in enumerate(attempts):
		# print(a)
		directory = dest_dir_pattern.replace("$a$",f"{attempt}")
		for file in glob.glob(directory+"297_*"):
			file_dir = "/".join(file.split('/')[:-1])
			file_name = file.split('/')[-1]
			# print(f"{file_dir}")
			os.system(f"mv {file} {file_dir}/{file_name.replace('297_','300_')}")
			# print(f"mv {file} {file_dir}/{file_name.replace("297","300")}")
			# print(file_name)
		# exit(0)




if __name__ == '__main__':
	main()