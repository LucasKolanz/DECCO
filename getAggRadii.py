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

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u



def calc_agg_radii(pos,COM):
	
	radii = 0
	for p in pos:
		r = np.linalg.norm(p-COM)
		
		if r > radii:
			radii = r 

	return radii



def main():

	#Open SpaceLab default file for directory information
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job_templates = [input_json["data_directory"] + 'jobsCosine/' + 'lognorm_relax' + '{a}/N_{n}/T_{t}/']
	job_templates.append(input_json["data_directory"] + 'jobsNovus/' + 'const_relax' + '{a}/N_{n}/T_{t}/')

	attempts = [i for i in range(30)]

	N = [30,100,300]
	# N=[300]

	Temps = [3,10,30,100,300,1000]
	# Temps = [1000]

	radii = np.full(shape=(len(job_templates)*len(N)*len(Temps)*len(attempts)),fill_value=np.nan)
	i = 0

	for job_template in job_templates:
		for a in attempts:
			for n in N:
				for t in Temps:
					folder = job_template.replace("{a}",str(a)).replace("{n}",str(n)).replace("{t}",str(t))
					
					try:
						pos,_,mass,_ = u.get_data(folder,relax=True)

						radii[i] = calc_agg_radii(pos,u.calcCOM(pos,mass))
						# print(radii[i])
					except:
						print(folder)

					i+=1
					
	print(np.nanmax(radii))
	print(np.nanmin(radii))

if __name__ == '__main__':
	main()