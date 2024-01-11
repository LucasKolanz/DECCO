import os
import glob
import numpy as np
import subprocess
import check_for_errors as cfe




def main():
	curr_folder = os.getcwd() + '/'

	job_folder = 'jobsCosine/'##FOR LOCAL
	# job_folder = 'jobs/'###FOR COSINE
	job_folder = 'erroredJobs/'

	job = curr_folder + job_folder + 'lognorm$a$/N_$n$/T_$t$/'
	# move_folder = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[30]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	error2_folders = cfe.check_final_error(job,cfe.ck_error2_by_file,N,Temps,attempts)

	print(error2_folders)


if __name__ == '__main__':
	main()