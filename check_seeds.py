import os
import glob
import numpy as np
import subprocess
import check_for_errors as cfe

def transfer(source,dest):
	if not os.path.exists(dest):
		os.makedirs(dest)
	os.system(f"cp {source}* {dest}.")

def get_attempt(folder):
	folder_base = folder 
	i = 1000
	folder = folder_base.replace("$a$",str(i))
	# print(folder+"timing.txt")
	# print(os.path.exists(folder+"timing.txt"))
	while os.path.exists(folder+"timing.txt"):
		i += 1
		folder = folder_base.replace("$a$",str(i))
		# print(folder+"timing.txt")
		# print(os.path.exists(folder+"timing.txt"))

	# print(folder)
	return i


def main():
	curr_folder = os.getcwd() + '/'

	job_folder = 'jobsCosine/'##FOR LOCAL

	base_job = curr_folder + job_folder + 'lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[30]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	for n in N:
		for t in Temps:
			print(f"====================  n={n} and t={t}  ====================")
			for a in attempts:
				folder = base_job.replace("$a$",str(a)).replace("$n$",str(n)).replace("$t$",str(t))
				
				if os.path.exists(folder+"seedFile.txt"):
					with open(folder+"seedFile.txt",'r') as f:
						print(f.read())


if __name__ == '__main__':
	main()