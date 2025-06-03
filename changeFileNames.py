"""
This file was originally written for SpaceLab/DECCO to check the seeds of all jobs in folders matching a specified pattern.

Author: Lucas Kolanz

This file is meant to go through all jobs contained in a specified file pattern and print out their seeds. All of these seeds should
be different, and this file is meant to help verify that.

TODO: use the seeds from input.json, not seedFile.txt. Save all seeds in a list to verify they are all unique.

"""




import h5py	
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



def rename_files(directory_in_str,just_print=True):

	#loop through files that for sure were from either the original naming convention
	#or start indexing at 0 rather than the number of balls in the sim

	#check if the zeroth index of either convention exists. If so, we need to rename files in this directory


	needsUpdate = False
	maxIndex = -1
	# path = '/'.join(directory_in_str.split('/')[0:-1]) + '/'

	if not os.path.isdir(directory_in_str):   
		print(f"Directory doesnt exist: {directory_in_str}")
		return

	directory = os.fsencode(directory_in_str)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.split("_")[0].isdigit():
			index = int(filename.split("_")[0])
			if index > maxIndex:
				maxIndex = index
		if filename.endswith("constants.csv") or filename.endswith("energy.csv") or filename.endswith("simData.csv"): 
			# if filename.startswith("2_R") or filename.startswith("0_") or not verify_index_based_on_balls(directory_in_str,index):
			if not verify_index_based_on_balls(directory_in_str,index):
				needsUpdate = True
			if filename.startswith("2_R"):
				file_base = "_".join(filename.split("_")[:-1])


	if needsUpdate:
		print("Updating directory: "+directory_in_str)
		# for i in [126,125,124,123,122]:
		for i in range(maxIndex,-1,-1): #This needs to go backwards so we don't overwrite anything. This is very important

			if i == 0:
				#newname is present for 0
				if os.path.exists(directory_in_str+"0_simData.csv"):
					glob_files = glob.glob(directory_in_str+"0_*")
				else:
					glob_files = glob.glob(directory_in_str+"2_R*")
			elif i == 2:
				#newname is present for 2
				if os.path.exists(directory_in_str+"2_simData.csv"):
					glob_files = glob.glob(directory_in_str+'2_*')
				else:
					glob_files = glob.glob(directory_in_str+'2_2_R*')

			else:
				glob_files = glob.glob(directory_in_str+str(i)+'_*')

			for full_path in glob_files:
				file = full_path.split("/")[-1]
				# prefix = "_".join(file.split("_")[:-1])
				suffix = file.split("_")[-1]
				newfile = str(i+3)+'_'+suffix

				newfull_path = directory_in_str + newfile

				#check we won't be overwriting any files with this name change
				# print('=======================')
				if not os.path.exists(newfull_path):
					if not verify_index_based_on_balls(directory_in_str,i):
						# pass
						if just_print:
							print("====================================================")
							print(f"original path: {full_path}")
							print(f"new path:      {newfull_path}")
						else:
							os.rename(full_path,newfull_path)

				else:
					print(f"WARNING already exists: {newfull_path}")




def rename_h5_files(directory_in_str,just_print=True):

	#loop through files that for sure were from either the original naming convention
	#or start indexing at 0 rather than the number of balls in the sim

	#check if the zeroth index of either convention exists. If so, we need to rename files in this directory


	needsUpdate = False
	maxIndex = -1
	# path = '/'.join(directory_in_str.split('/')[0:-1]) + '/'

	if not os.path.isdir(directory_in_str):   
		print(f"Directory doesnt exist: {directory_in_str}")
		return

	directory = os.fsencode(directory_in_str)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.split("_")[0].isdigit():
			index = int(filename.split("_")[0])
			if index > maxIndex:
				maxIndex = index
		# if filename.endswith(".h5"): 
		# 	# if filename.startswith("2_R") or filename.startswith("0_") or not verify_index_based_on_balls(directory_in_str,index):
		# 	if not verify_index_based_on_balls(directory_in_str,index):
		# 		needsUpdate = True
		# 	# if filename.startswith("2_R"):
			# 	file_base = "_".join(filename.split("_")[:-1])


	# if needsUpdate:
	# for i in [126,125,124,123,122]:
	for i in range(maxIndex,-1,-1): #This needs to go backwards so we don't overwrite anything. This is very important

		# if i == 0:
		# 	#newname is present for 0
		# 	if os.path.exists(directory_in_str+"0_simData.csv"):
		# 		glob_files = glob.glob(directory_in_str+"0_*")
		# 	else:
		# 		glob_files = glob.glob(directory_in_str+"2_R*")
		# elif i == 2:
		# 	#newname is present for 2
		# 	if os.path.exists(directory_in_str+"2_simData.csv"):
		# 		glob_files = glob.glob(directory_in_str+'2_*')
		# 	else:
		# 		glob_files = glob.glob(directory_in_str+'2_2_R*')

		# else:
		glob_files = glob.glob(directory_in_str+str(i)+'_*')

		for full_path in glob_files:
			file = full_path.split("/")[-1]
			# prefix = "_".join(file.split("_")[:-1])
			suffix = file.split("_")[-1]
			newfile = str(i+3)+'_'+suffix

			newfull_path = directory_in_str + newfile

			#check we won't be overwriting any files with this name change
			# print('=======================')
			if not verify_index_based_on_balls(directory_in_str,i):
				if not os.path.exists(newfull_path):
						# pass
						if just_print:
							print("====================================================")
							print(f"original path: {full_path}")
							print(f"new path:      {newfull_path}")
						else:
							os.rename(full_path,newfull_path)

				else:
					print(f"WARNING already exists: {newfull_path}")
	# exit(0)


				
def verify_directory_based_on_balls(directory_in_str):

	#loop through files that for sure were from either the original naming convention
	#or start indexing at 0 rather than the number of balls in the sim

	#check if the zeroth index of either convention exists. If so, we need to rename files in this directory


	needsUpdate = False
	maxIndex = -1
	minIndex = float("inf")
	# path = '/'.join(directory_in_str.split('/')[0:-1]) + '/'

	if os.path.isdir(directory_in_str):    
		directory = os.fsencode(directory_in_str)
		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			# print(filename)
			if filename.split("_")[0].isdigit():
				index = int(filename.split("_")[0])
				if index > maxIndex:
					maxIndex = index
				if index < minIndex:
					minIndex = index
			# if filename.endswith("constants.csv") or filename.endswith("energy.csv") or filename.endswith("simData.csv"): 
			# 	if filename.startswith("2_R") or filename.startswith("0_"):
			# 		needsUpdate = True
			# 	if filename.startswith("2_R"):
			# 		file_base = "_".join(filename.split("_")[:-1])

		

		print("checking directory: "+directory_in_str)
		# for i in [126,125,124,123,122]:
		if maxIndex != -1 and minIndex != float("inf"):
			for i in range(maxIndex,minIndex-1,-1): #This needs to go backwards so we don't overwrite anything. This is very important
				#check the constants file
				constants_file = directory_in_str+str(i)+"_constants.csv"
				if os.path.exists(constants_file):
					if len(u.get_radii(directory_in_str,i,relax=False)) != i:
						print(f"index doesn't match ball number for index {i}")
				# else:
				# 	print(f"WARNING:: FILE MISSING for index {i}")

				#check the simData file
				simData_file = directory_in_str+str(i)+"_simData.csv"
				if os.path.exists(simData_file):
					input_file = directory_in_str+'input.json'
					with open(input_file,'r') as fp:
						input_json = json.load(fp)
					simTimeSeconds = input_json['simTimeSeconds']
					timeResolution = input_json['timeResolution']
					numLines = int(simTimeSeconds/timeResolution)

					with open(simData_file, 'r') as file:
						length = len(file.readlines())
					if length > numLines+3:
						print(f"simData file {simData_file} is {length} lines. This is bad.")
				# else:
				# 	print(f"WARNING:: FILE MISSING for index {i}")

def get_num_balls(directory_in_str,index):
	relax = False
	if "relax" in directory_in_str:
		relax = True
	return len(u.get_radii(directory_in_str,index,relax=relax))

def verify_index_based_on_balls(directory_in_str,index):
	relax = False
	if "relax" in directory_in_str:
		relax = True

	if get_num_balls(directory_in_str,index) != index:
		return False
	else:
		return True

		




def main():

	#Open SpaceLab default file for directory information
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_templates = [input_json["data_directory"] + 'jobsCosine/' + 'lognormrelax_{a}/N_{n}/T_{t}/']
	job_templates = [input_json["data_directory"] + 'jobsNovus/' + 'const_{a}/N_{n}/T_{t}/']
	job_templates = [input_json["data_directory"] + 'jobs/' + 'BAPA_{a}/M_{m}/N_{n}/T_{t}/']


	attempts = [i for i in range(10)]
	# attempts = [i for i in range(30,40)]

	# attempts = [0]

	M = [3,5,10,15]

	# N = [30,100,300]
	N=[300]

	# Temps = [3,10,30,100,300,1000]
	Temps = [1000]

	# i = 0

	just_print = True

	for job_template in job_templates:
		for a in attempts:
			for m in M:
				for n in N:
					for t in Temps:
						folder = job_template.replace("{a}",str(a)).replace("{m}",str(m)).replace("{n}",str(n)).replace("{t}",str(t))
						# print(folder)
						# rename_files(folder,just_print)
						# rename_h5_files(folder,just_print)
						verify_directory_based_on_balls(folder)


			# for n in N:
			# 	for t in Temps:
			# 		folder = job_template.replace("{a}",str(a)).replace("{n}",str(n)).replace("{t}",str(t))
			# 		print(folder)
			# 		rename_files(folder,just_print)
			# 		# rename_h5_files(folder,just_print)
					
	# if not just_print:
	# 	for job_template in job_templates:
	# 		for a in attempts:
	# 			for n in N:
	# 				for t in Temps:
	# 					folder = job_template.replace("{a}",str(a)).replace("{n}",str(n)).replace("{t}",str(t))
	# 					verify_name_based_on_balls(folder)

if __name__ == '__main__':
	main()