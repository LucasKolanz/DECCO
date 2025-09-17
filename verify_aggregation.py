"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file goes through all folders in a specified job set and verifys the state of all indices.
It gets the final state of each index and ensures that all balls are included in the final aggregate

"""





import sys
import glob
import os
import json
import numpy as np
import argparse

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u


def are_spheres_connected(pos, radii):
	from collections import deque

	def distance(a, b):
		return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5

	def is_connected(index1, index2):
		return distance(pos[index1], pos[index2]) <= 1.1*radii[index1] + 1.1*radii[index2]

	n = len(pos)
	if n == 0:
		return False

	# Build the adjacency list
	adj_list = {i: [] for i in range(n)}
	for i in range(n):
		for j in range(i + 1, n):
			if is_connected(i, j):
			    adj_list[i].append(j)
			    adj_list[j].append(i)

	# BFS to check connectivity
	visited = [False] * n
	queue = deque([0])
	visited[0] = True
	while queue:
		node = queue.popleft()
		for neighbor in adj_list[node]:
			if not visited[neighbor]:
				visited[neighbor] = True
				queue.append(neighbor)

	# print(visited)
	return all(visited)


def delete_if_exists(file,actually_delete=False):
	if os.path.exists(file):
		if actually_delete:
			os.remove(file)
		else:
			print(f"file to delete: {file}")

# def stop_job(directory):



if __name__ == '__main__':
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobs/BAPA_*'
	# data_prefolder = path + 'jobs/constrollingfric*'
	# data_prefolder = '/home/kolanzl/Desktop/Visualize/V19/'

	possible_dirs = u.get_directores(data_prefolder)

	failed_folders = []
	failed_files = []
	attErrs = []


	for directory in possible_dirs:
		
		has_failed = False
		failed_indices = []
		for index in u.get_all_indices(directory,checkpoint=False):
			# print(index)
			try:
				pos,_,_,radius,_,_ = u.get_all_data(directory,data_index=index,linenum=-1,relax=False)
			except:
				pos = []
				radius = []
				attErrs.append(directory)
			if len(pos) > 0 and not are_spheres_connected(pos,radius):
				has_failed = True
				failed_indices.append(index)
		if has_failed:
			print(f"directory failed connection check: {directory}")
			print(f"\tfor indices: {failed_indices}")
			# failed_files.append([str(i)+"_*" for i in failed_indices])
			# failed_folders.append(directory)
	
	print("The following folders had AttributeError when attempting to get data:")
	print(attErrs)

	print(failed_folders)
	print(failed_files)

	# for f_i,failed_folder in enumerate(failed_folders):
	# 	print(f"folder: {failed_folder}")
	# 	for file in failed_folder[f_i]:
	# 		print(f"\t{file}")


	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--delete', default=False, help="True to delete failed runs, False to list them")
	# parser.add_argument('-n', '--novus', default=False, help="True to check if job is running and to stop it if it is. False to skip running check")

	args = parser.parse_args()

	actually_delete = args.delete

	print("The following files contain failed aggregation:")
	for f_i,failed_folder in enumerate(failed_folders):
		print(f"failed folder: {failed_folder}")

		# if args.novus:
		# 	stop_job(failed_folder)

		delete_if_exists(failed_folder+"timing.txt",actually_delete)
		delete_if_exists(failed_folder+"job_data.csv",actually_delete)

		for file in failed_files[f_i]:
			for globfile in glob.glob(failed_folder+file):
				delete_if_exists(globfile,actually_delete)
