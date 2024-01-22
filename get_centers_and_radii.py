import os
import glob
import numpy as np
import sys
import json
import h5py

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u


def get_last_line(folder_path,data_index=-1):
	file_path = u.get_data_file(folder_path,data_index)
	with open(folder_path+file_path, 'rb') as f:
		try:  # catch OSError in case of a one line file 
			f.seek(-2, os.SEEK_END)
			while f.read(1) != b'\n':
				f.seek(-2, os.SEEK_CUR)
		except OSError:
			f.seek(0)
		last_line = f.readline().decode()
	return last_line

def get_last_line_data(data_folder,data_index=-1): #Works with csv and h5
	# data_headers = np.loadtxt(data_folder + data_file,skiprows=0,dtype=str,delimiter=',')[0]
	csv_data = False
	h5_data = False
	data_file = u.get_data_file(data_folder,data_index)
	if data_file.endswith(".csv"):
		csv_data = True
	elif data_file.endswith(".h5"):
		h5_data = True
	# print("data file: {}".format(data_file))
	if csv_data:
		try:
			data = np.loadtxt(data_folder + data_file,skiprows=1,dtype=float,delimiter=',')
			if data.ndim > 1:
				data = data[-1]
			# print(data)
			# print(data_folder + data_file)
		except Exception as e:
			with open(data_folder + data_file) as f:
				for line in f:
					pass
				last_line = line
			data = np.array([last_line.split(',')],dtype=np.float64)
			# print(data)
			print("ERROR CAUGHT getting data in folder: {}".format(data_folder))
			print(e)
			return None
	elif h5_data:
		data = u.get_last_line_h5data_from_file(data_folder+data_file)
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")
	# print("DATA LEN: {} for file {}{}".format(data.size,data_folder,data_file))
	# print("FOR {} Balls".format(data.size/11))
	return data

def get_constants(data_folder,data_index=-1):
	csv_data = False
	h5_data = False

	data_file = u.get_data_file(data_folder,data_index)

	if data_file.endswith(".csv"):
		csv_data = True
		data_file = data_file.replace("simData","constants")
	elif data_file.endswith(".h5"):
		h5_data = True
	# print("data file: {}".format(data_file))
	if csv_data:
		data = np.loadtxt(data_folder + data_file,skiprows=0,dtype=float,delimiter=',').reshape(-1)
		# try:
		# 	# print(data)
		# 	# print(data_folder + data_file)
		# except Exception as e:
		# 	with open(data_folder + data_file) as f:
		# 		for line in f:
		# 			pass
		# 		last_line = line
		# 	data = np.array([last_line.split(',')],dtype=np.float64)
		# 	# print(data)
		# 	print("ERROR CAUGHT getting data in folder: {}".format(data_folder))
		# 	print(e)
		# 	return None
	elif h5_data:
		# data = u.get_last_line_h5data_from_file(data_folder+data_file)
		with h5py.File(data_folder+data_file, 'r') as file:
			data = file['/constants'][:]
			# data = np.array(consts).reshape(int(len(consts)/3),3)
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")
	# print("DATA LEN: {} for file {}{}".format(data.size,data_folder,data_file))
	# print("FOR {} Balls".format(data.size/11))
	return data

def center_radii(job,n):

	simData = get_last_line_data(job,n-1)


	const_data = get_constants(job,n-1)
	# const_data = np.loadtxt(base_file_path+'constants.csv',delimiter=',',dtype=np.float64)
	radii = np.array([const_data[i] for i in range(len(const_data)) if i % 3==0])
	# print(radii)

	output = []

	for i in range(0,len(simData),11):
		output.append([simData[i],simData[i+1],simData[i+2],radii[int(i/11)]])

	return np.array(output)


def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	# job = curr_folder + 'jobs/' + job_set_name + str(attempt) + '/'
	data_dir = input_json["data_directory"]


	job_set_name = "const"
	job_name = "jobsNovus"

	job_set_name = "lognorm"
	job_name = "jobsCosine"

	attempts = [i for i in range(30)]
	# attempts = [1]

	# N = [30,100]
	N=[30]
	# N=[10]

	# Temps = [3,10,30,100,300,1000]
	Temps = [3]


	for n in N:
		for t_i,Temp in enumerate(Temps):
			output = np.full(shape=(len(attempts),n+2,4),fill_value=np.nan,dtype=np.float64)
			output_file_name = f"{data_dir}data/center_radii_{job_name}_N-{n}_T-{Temp}.csv"
			for a_i,attempt in enumerate(attempts):
				job_folder = data_dir + job_name + '/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				# job = job_folder+str(n-1)+"_2_*"
				
				if os.path.exists(job_folder):
					cr_output = center_radii(job_folder,n)
					if not isinstance(cr_output,int):
						output[a_i,:,:] = cr_output[:,:]

			np.savetxt(output_file_name,output.reshape(-1))

			#test output
			check_ouput = np.loadtxt(output_file_name,dtype=np.float64).reshape(len(attempts),n+2,4)

			if (np.array_equal(check_ouput,output,equal_nan=True)):
				print("SUCCESS")
				print(check_ouput.shape)
			else:
				print("FAILURE")
	# print(output)

if __name__ == '__main__':
	main()

#/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/29_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k3e+00_Ha5e-12_dt4e-10_simData.csv
