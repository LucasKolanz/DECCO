##This is a file with helpful python functions for this project

import os
import json
import numpy as np
import h5py	

def index_from_file(file):
	file_split = file.split("_")
	if not file_split[1].isnumeric() and len(file_split)>2:
		return 0
	else:
		return int(file_split[0])

def find_max_index(folder,relax=False):
	max_index = -1
	rel = ""
	if relax:
		rel = "RELAX"

	files = os.listdir(folder)
	for file in files:
		if file.split("_")[0].isnumeric():
			if file.endswith(f"{rel}simData.csv") or file.endswith(f"{rel}constants.csv") or file.endswith(f"{rel}energy.csv") or file.endswith("data.h5"):
				max_index = max(index_from_file(file),max_index)
	return max_index

def getSimDataFile(folder,relax=False,index=-1):
	return getDataFile(folder,relax,index)

def getConstantsFile(folder,relax=False,index=-1):
	fileBase = getDataFile(folder,relax,index)
	if fileBase.endswith(".h5"):
		return fileBase
	elif fileBase.endswith(".csv"):
		return fileBase.replace("simData","constants")
	else:
		return None

def getEnergyFile(folder,relax=False,index=-1):
	fileBase = getDataFile(folder,relax,index)
	if fileBase.endswith(".h5"):
		return fileBase
	elif fileBase.endswith(".csv"):
		return fileBase.replace("simData","energy")
	else:
		return None


def getDataFile(folder,relax=False,index=-1):
	rel = ""
	if relax:
		rel = "RELAX"
	if (os.path.exists(folder+"timing.txt")):
		directory = os.fsencode(folder)

		if index < 0:
			index = find_max_index(folder,relax)
			if index < 0:
				print(f"ERROR: find_max_index found no indicies in folder {folder}")


		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			if filename.startswith(str(index)) and (filename.endswith(f"{rel}simData.csv") or filename.endswith(f"{rel}data.h5")):
				return filename
	else:
		print(f"ERROR: File does not exist: {folder}timing.txt")

	print(f"ERROR: No {rel}simData file in {folder} with index {index}")

	return None


def getProperties(global_path):
	with open("/".join(global_path.split("/")[:-1])+"/input.json",'r') as fp:
		input_json = json.load(fp)
		return input_json["properties"]

	print(f"ERROR: Couldn't find properties for {global_path}")
	return -1

def getConstants(global_path):
	if global_path.endswith(".csv"):
		data_constants = np.loadtxt(global_path,skiprows=0,dtype=float,delimiter=',')
	elif global_path.endswith(".h5"):
		with h5py.File(global_path, 'r') as file:
			consts = file['/constants'][:]
			data_constants = np.array(consts).reshape(int(len(consts)/3),3)
	else:
		print(f"ERROR: data file type not recognized by utils.py: {data_file}")
	return data_constants

def getSimData(global_path):
	if global_path.endswith(".csv"):
		return getPosCSV(global_path)
	elif global_path.endswith(".h5"):
		return getPosH5(global_path)
	else:
		print(f"ERROR: file not of type csv or h5 '{global_path}'")
	return simData

def getPosCSV(global_path):
	simData_properties = getProperties(global_path)
	simData = np.loadtxt(global_path,skiprows=1,delimiter=',',dtype=np.float64)[-1]
	simData = simData.reshape(int(simData.shape[0]/simData_properties),simData_properties) #now it is simData[ball,property]
	return simData[:,:3]

def getPosH5(global_path):
	with h5py.File(global_path, 'r') as file:
		constants = np.array(file['/constants'][:])
		radii = constants[np.where(np.arange(constants.shape[0])%3==0)]
		num_spheres = radii.shape[0]

		simData_single_ball_width = getProperties(global_path)

		simData = np.array(file['/simData'][-simData_single_ball_width*num_spheres:])
		simData = simData.reshape(num_spheres,simData_single_ball_width)
		pos = simData[:,:3]
	return pos

def getPos(folder,relax=False,index=-1):
	file = getSimDataFile(folder,relax,index)
	if file is None:
		return None
	simData = getSimData(folder+file)
	return simData[:,:3]


def getMass(folder,relax=False,index=-1):
	file = getConstantsFile(folder,relax,index)
	if file is None:
		return None
	consts = getConstants(folder+file)
	return consts[:,0]


def verifyExtractPos(call, result):
	if result is None:
		print(f'FAILED: {call}')
	else:
		print(f'PASSED: {call}')

if __name__ == '__main__':
	pos = getPos("/media/kolanzl/easystore/SpaceLab_data/jobs/relaxTest6/N_300/T_3/")
	verifyExtractPos('getPos("/media/kolanzl/easystore/SpaceLab_data/jobs/relaxTest6/N_300/T_3/")',pos)

	pos = getPos("/media/kolanzl/easystore/SpaceLab_data/jobs/relaxTest6/N_300/T_3/",True,297)
	verifyExtractPos('getPos("/media/kolanzl/easystore/SpaceLab_data/jobs/relaxTest6/N_300/T_3/"),True,297',pos)
	
	pos = getPos("/media/kolanzl/easystore/SpaceLab_data/jobsNovus/const0/N_300/T_3/")
	verifyExtractPos('getPos("/media/kolanzl/easystore/SpaceLab_data/jobsNovus/const0/N_300/T_3/")',pos)

	pos = getPos("/media/kolanzl/easystore/SpaceLab_data/jobsNovus/const0/N_300/T_3/",False,60)
	verifyExtractPos('getPos("/media/kolanzl/easystore/SpaceLab_data/jobsNovus/const0/N_300/T_3/"),False,60',pos)
	