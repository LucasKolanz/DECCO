import os
import json
import multiprocessing as mp
import subprocess

def run_job(location,num_balls):
	cmd = ["python3", "{}run_sim.py".format(location), location, str(num_balls)]
	print('run cmd', cmd)
	subprocess.run(cmd)

if __name__ == '__main__':
	#make new output folders
	curr_folder = os.getcwd() + '/'

	try:
		# os.chdir("{}ColliderSingleCore".format(curr_folder))
		subprocess.run(["make","-C","ColliderSingleCore"], check=True)
	except:
		print('compilation failed')
		exit(-1)

	job_set_name = "vary_g_size"
	# folder_name_scheme = "g_"

	#Make an array of what you want to vary

	runs_at_once = 2
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 
	attempts = [1] 
	grid_sizes = [2e-5,6e-5,1e-4,1.4e-4,1.8e-4,2e-4]
	N = 100
	Temps = [3]
	folders = []
	for attempt in attempts:
		for gs in grid_sizes:
			for Temp in Temps:
				job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'gs_' + str(gs) + '/'
				if not os.path.exists(job):
					os.makedirs(job)
				else:
					print("Job '{}' already exists.".format(job))


				#load default input file
				with open(curr_folder+"default_files/default_input.json",'r') as fp:
					input_json = json.load(fp)

				####################################
				######Change input values here######
				input_json['temp'] = Temp
				input_json['gridSize'] = gs
				####################################

				with open(job + "input.json",'w') as fp:
					json.dump(input_json,fp,indent=4)

				###Copy over some data so we don't start from beginning
				cp_source_folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/tempVarianceRand_attempt1/T_3/'
				cp_source_1 = '80_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_'
				cp_source_2 = '81_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_'
				for source in [cp_source_1,cp_source_2]:
					for file in ['simData.csv','constants.csv','energy.csv']:
						os.system("cp {}{}{} {}{}{}".format(cp_source_folder,source,file,job,source,file))

				#add run script and executable to folders
				os.system("cp default_files/run_sim.py {}run_sim.py".format(job))
				os.system("cp ColliderSingleCore/ColliderSingleCore.o {}ColliderSingleCore.o".format(job))
				folders.append(job)
	# print(folders)
	# if len(N) != len(folders):
	N = [str(N) for i in range(len(folders))]

	inputs = list(zip(folders,N))
	# print(inputs)

	for i in range(0,len(folders),runs_at_once):
		with mp.Pool(processes=runs_at_once) as pool:
			pool.starmap(run_job,inputs[i:i+runs_at_once]) 


	
