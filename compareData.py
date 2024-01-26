import os
import numpy as np
import sys
import h5py
relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'
sys.path.append(project_path+"utilities/")
import utils as u
# import porosity_FD as p

def check_arrays(array1,array2,tol):
	different = np.zeros(array1.shape,dtype=bool)
	for i in range(len(array1)):
		if np.average([array1[i],array2[i]]) != 0:
			if np.abs(array1[i]-array2[i])/np.average([array1[i],array2[i]]) > tol:
				different[i] = 1
			
	return different

def get_data(file,dtype=""):
	if file.endswith(".csv"):
		if dtype == "simData" or dtype == "energy":
			skip = 1
		elif dtype == "constants":
			skip = 0

		file = file.replace("simData",dtype)

		# print(np.loadtxt(file,delimiter=",",dtype=np.float64,skiprows=skip).reshape(-1).shape)
		return np.loadtxt(file,delimiter=",",dtype=np.float64,skiprows=skip).reshape(-1)
	elif file.endswith(".h5"):
		with h5py.File(file,'r') as f:
			return np.array(f[f'/{dtype}'][:])

# def get_old_csv()

def main():
	folder1 = "/home/lucas/Desktop/SpaceLab_data/oldVersionCSVData/"
	folder2 = "/home/lucas/Desktop/SpaceLab_data/jobs/errorckcsv1/N_10/T_3/"
	folder1 = "/global/homes/l/lpkolanz/DECCO/SpaceLab_data/jobs/masterComp5/N_10/T_3/" #Compare to 5 because new calculation order
	folder2 = "/global/homes/l/lpkolanz/DECCO/SpaceLab_data/jobs/masterComp13/N_10/T_3/"

	N = 10
	temp = 100
	body = ""
	old1 = False;
	old2 = False;
	body = "2_R6e-05_v4e-01_cor0.63_mu0.1_rho2.25_k5e+00_Ha5e-12_dt4e-10_"
	tol = 1e-5
	for ind in range(N):

		file1 = u.get_data_file(folder1,ind,old1)
		file2 = u.get_data_file(folder2,ind,old2)
		# print(folder1+file1)
		# print(folder2+file2)

		print(f"===================TESTING simData {ind}===================")
		simData1 = get_data(folder1+file1,"simData")#np.array(f['/simData'][:])
		simData2 = get_data(folder2+file2,"simData")#np.loadtxt(folder2+file2+"simData.csv",delimiter=",",dtype=np.float64,skiprows=1).reshape(simData1.shape)
		different = check_arrays(simData1,simData2,tol);
		# different = ~np.isclose(simData1,simData2,tol);
		if np.sum(different) == 0:
			print("No errors")
		else:
			print(different[:20])
			print(simData1[:20])
			print(simData2[:20])
			print(f"{np.sum(different)}/{simData1.shape} different values")
		print(f"===================TEST Finished===================")
		print(f"===================TESTING constants {ind}===================")
		constData1 = get_data(folder1+file1,"constants")#np.array(f['/constants'][:])
		constData2 = get_data(folder2+file2,"constants")#np.loadtxt(folder2+file2+"constants.csv",delimiter=",",dtype=np.float64,skiprows=0).reshape(constData1.shape)
		different = check_arrays(constData1,constData2,tol);
		# different = ~np.isclose(constData1,constData2,tol);
		if np.sum(different) == 0:
			print("No errors")
		else:
			print(constData1[different])
			print(constData2[different])
			print(f"{np.sum(different)}/{constData1.shape} different values")
		print(f"===================TEST Finished===================")
		print(f"===================TESTING energy {ind}===================")
		energyData1 = get_data(folder1+file1,"energy")#np.array(f['/energy'][:])
		energyData2 = get_data(folder2+file2,"energy")#np.loadtxt(folder2+file2+"energy.csv",delimiter=",",dtype=np.float64,skiprows=1).reshape(energyData1.shape)
		different = check_arrays(energyData1,energyData2,tol);
		# different = ~np.isclose(energyData1,energyData2,tol);
		if np.sum(different) == 0:
			print("No errors")
		else:
			# print(energyData1[different])
			# print(energyData2[different])
			print(f"{np.sum(different)}/{energyData1.shape} different values")
		print(f"===================TEST Finished===================")

if __name__ == '__main__':
	main()
