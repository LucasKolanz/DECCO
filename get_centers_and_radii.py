import os
import glob
import numpy as np


def get_last_line(file_path):
	with open(file_path, 'rb') as f:
		try:  # catch OSError in case of a one line file 
			f.seek(-2, os.SEEK_END)
			while f.read(1) != b'\n':
				f.seek(-2, os.SEEK_CUR)
		except OSError:
			f.seek(0)
		last_line = f.readline().decode()
	return last_line

def center_radii(job,n):
	file_paths = glob.glob(job)
	base_file_path = '_'.join(file_paths[0].split('_')[:-1])+'_'

	last_line = get_last_line(base_file_path+'simData.csv')
	const_data = np.loadtxt(base_file_path+'constants.csv',delimiter=',',dtype=np.float64)
	radii = const_data[:,0]

	output = []
	
	line_split = last_line.split(',')

	if len(line_split)/11 == 1.0*(n+2):
		for i in range(0,len(line_split),11):
			output.append([line_split[i],line_split[i+1],line_split[i+2],radii[int(i/11)]])
		return np.array(output)
	else:
		print(f"data corrupt in: {job}")
		return -1


def main():

	curr_folder = os.getcwd() + '/'

	job_set_name = "lognorm"

	attempts = [i for i in range(30)]

	N = [30,100]
	# N=[100]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]


	for n in N:
		for t_i,Temp in enumerate(Temps):
			output = np.full(shape=(len(attempts),n+2,4),fill_value=np.nan,dtype=np.float64)
			output_file_name = f"{curr_folder}data/center_radii_N-{n}_T-{Temp}.csv"
			for a_i,attempt in enumerate(attempts):
				job_folder = curr_folder + 'jobsCosine/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				job = job_folder+str(n-1)+"_2_*"
				
				cr_output = center_radii(job,n)
				# print(cr_output.shape)
				if not isinstance(cr_output,int):
					output[a_i,:,:] = cr_output[:,:]

			np.savetxt(output_file_name,output.reshape(-1))

			#test output
			check_ouput = np.loadtxt(output_file_name,dtype=np.float64).reshape(len(attempts),n+2,4)
			if (np.array_equal(check_ouput,output)):
				print("SUCCESS")
				print(check_ouput.shape)

if __name__ == '__main__':
	main()

#/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/29_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k3e+00_Ha5e-12_dt4e-10_simData.csv
