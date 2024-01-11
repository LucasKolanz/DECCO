from __future__ import division
import bpy
import numpy as np
from mathutils import *
from math import *
import fnmatch
import os
import h5py

def get_simData_and_consts(path,filename,fileindex,csv):
	numSpheres = -1
	steps = -1
	
	if csv:
		filename = "_"+"_".join(filename.split('_')[1:-1]) +'_'
		print(filename)
		
		print(path + str(fileindex) + filename + "simData.csv")
		if fileindex == 0:
			if "SpaceLab_branch" in path.split("/"):
				pth = path + filename[2:] + "simData.csv"
			else:
				pth = path + filename[1:] + "simData.csv"
		else:
			if "SpaceLab_branch" in path.split("/"):
				pth = path + str(fileindex) + filename[:-1] + "simData.csv"
			else:
				pth = path + str(fileindex) + filename + "simData.csv"
		try:
			simData = np.loadtxt(pth,dtype=float,delimiter=',',skiprows = 1)
		except Exception as e:
			print("ERROR CAUGHT in folder: {}".format(pth))
			print(e)
		#    simData = np.array([last_line.split(',')],dtype=np.float64)
		print("PATH")
		print(pth)
		print("DATA")
		print(simData)
		#simData = np.array([simData]) # Uncomment this line for single timestep data with no headers
		#simData.T # Uncomment this line for single timestep data with no headers
		steps = len(simData)
		print("steps: ",steps)
		if fileindex == 0:
			if "SpaceLab_branch" in path.split("/"):
				constants = np.genfromtxt(path + filename[2:] + "constants.csv",dtype=float,delimiter=',')
			else:
				constants = np.genfromtxt(path + filename[1:] + "constants.csv",dtype=float,delimiter=',')
		else:
			if "SpaceLab_branch" in path.split("/"):
				constants = np.genfromtxt(path + str(fileindex) + filename[1:] + "constants.csv",dtype=float,delimiter=',')
			else:
				print("==================")
				print(path + str(fileindex) + filename + "constants.csv")
				print("==================")
				constants = np.genfromtxt(path + str(fileindex) + filename + "constants.csv",dtype=float,delimiter=',')
		numSpheres = len(constants)
	else:
		filename = '_'+filename.split('_')[-1]
		
		print(path+str(fileindex)+filename)
		
		with h5py.File(path+str(fileindex)+filename, 'r') as file:
			simData = np.array(file['simData'])
			constants = np.array(file['constants'])
			numSpheres = (int)(constants.size/3) #3 is the width of the constants file (see DECCOData)
			simRows = (int)(simData.size/(numSpheres*11)) #11 is the single ball width of simData (see DECCOData)
			steps = simRows
			simCols = numSpheres*11
			simData = simData.reshape(simRows,simCols)
			constants = constants.reshape(numSpheres,(int)(constants.size/numSpheres))
	return [simData,constants,numSpheres,steps]
	

# Delete old stuff first:
foo_objs = [obj for obj in bpy.context.scene.objects if fnmatch.fnmatchcase(obj.name, "*phere*")]

for obj in foo_objs:
	bpy.data.objects.remove(obj, do_unlink = True)

foo_objs = [obj for obj in bpy.context.scene.objects if fnmatch.fnmatchcase(obj.name, "Mball*")]

for obj in foo_objs:
	bpy.data.objects.remove(obj, do_unlink = True)
	

stepSkip = 1
stepTime = 1e-5
properties = 11
scaleUp = 1e5
frameNum = 0

#path = "C:/Users/milin/Desktop/VSOUT/1000-0/201_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/3-0/201_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/10-1e-4/228_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"

#path = "C:/Users/milin/Desktop/VSOUT/300/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/100/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/30/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/10/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"

#path = "C:/Users/milin/Desktop/VSOUT/1000-1e-4/50_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/300/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/100/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/30/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/10/213_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/3-1e-4/50_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
	 
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregateNN_30/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregateNN_30_1/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/tempVariance_attempt4/T_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregateNN_30_3/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregategridNN/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregategridNN_1/N_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/large_aggregate_optO3_1/N_1000/T_3/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/PairParallelTest1/N_1000/T_500/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/test1/N_10/T_1000/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error_checking1/N_1/T_1/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/cuttoff_test/c_1.2/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/step_back18/N_2/T_1/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/mu_max9/N_30/T_3/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/testingTime3/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/calibrateTest1/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/mu_max8/N_100/T_30/'
path = '/home/lpkolanz/Desktop/SpaceLab/jobs/accuracyTest1/'
path = '/home/lpkolanz/Desktop/SpaceLab/jobs/parTest/'
path = '/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab/jobs/accuracyTest5/N_10/T_100/'
path = '/home/lpkolanz/Desktop/SpaceLab/jobs/singleCoreComp/'
path = '/home/lpkolanz/Desktop/SpaceLab/jobs/singleCoreComparison_COPY7/'
path = '/home/lpkolanz/Desktop/SpaceLab/jobs/singleCoreComparison6/'
path = '/home/lpkolanz/Desktop/SpaceLab_branch_copy/SpaceLab/testing/jobs/multiCoreTest1/'
path = '/home/kolanzl/Desktop/SpaceLab_copy/SpaceLab_data/test1/N_5/T_3/'
path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/tempVarianceRand_attempt8/N_30/T_3/'
#filename = '_2_R5e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_'
#filename = '_1_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_' 
#filename = '_'
#filename = '_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-12_'
#filename = '_2_R2e+05_v4e-01_cor0.63_mu0.1_rho2.25_k2e+10_Ha5e-12_dt3e+00_'
#filename = '_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt4e-10_'
simStart = 0
simEnd = 2

csv = False
filename = ''
for file in os.listdir(path):
#    print(file)
	if len(file) > len(filename):
		if file[-4:] == ".csv":
			filename = file
			csv = True
		elif file[-3:] == ".h5":
			filename = file	
	
#print("num constants: ",numSpheres)

print(filename)

simData,constants,numSpheres,steps = get_simData_and_consts(path,filename,simStart,csv)

sphereSet = []
actionSet = []

# Initial sphere mesh to be instanced:
bpy.ops.mesh.primitive_uv_sphere_add(location = (0,0,0), radius = 1)
#bpy.ops.object.metaball_add(type='BALL', radius=2, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
obj = bpy.context.object # the currently selected object
#obj.data.resolution = .1
#obj.data.threshold = 1.7
sphereMesh = obj.data # retrieve the mesh
bpy.data.objects.remove(obj) # remove the object


# Instanciate spheres:
for sphere in range(numSpheres):
	sphereSet.append(bpy.data.objects.new("Mball." + str(sphere),sphereMesh))
	bpy.context.scene.collection.objects.link(sphereSet[sphere]) # link the object to the scene collection
	sphereSet[sphere].scale = (scaleUp*constants[sphere,0],scaleUp*constants[sphere,0],scaleUp*constants[sphere,0])
	sphereSet[sphere].location = (scaleUp*simData[0][0 + properties*sphere],scaleUp*simData[0][1 + properties*sphere],scaleUp*simData[0][2 + properties*sphere])
	sphereSet[sphere].rotation_mode = "XYZ"
#print(sphereSet)
for sim in range(simStart,simEnd+1):
	if frameNum > 1:
		print(simData.shape)        
		simData,constants,numSpheres,steps = get_simData_and_consts(path,filename,sim,csv)
		print(sim)
		print(simData.shape)
		print(numSpheres)
		
		# Instanciaten the new particle and the end of file:
		sphereSet.append(bpy.data.objects.new("Mball." + str(numSpheres),sphereMesh))
		print(len(sphereSet))
		sphere += 1
		bpy.context.scene.collection.objects.link(sphereSet[numSpheres-1]) # link the object to the scene collection
		sphereSet[sphere].scale = (scaleUp*constants[numSpheres-1,0],scaleUp*constants[numSpheres-1,0],scaleUp*constants[numSpheres-1,0])
		sphereSet[sphere].location = (scaleUp*simData[0][0 + properties*(numSpheres-1)],scaleUp*simData[0][1 + properties*(numSpheres-1)],scaleUp*simData[0][2 + properties*(numSpheres-1)]) 
		sphereSet[sphere].rotation_mode = "XYZ"

	for step in range(steps):
#        print("step: ",step) 
		if step % steps / 10 == 0:
			print(str(step/steps*100)+'%')
		if step % stepSkip == 0:
			bpy.context.scene.frame_set(frameNum)
			for sphere in range(numSpheres):
				
				# Move spheres
				sphereSet[sphere].location = (scaleUp*simData[step][0 + properties*sphere],scaleUp*simData[step][1 + properties*sphere],scaleUp*simData[step][2 + properties*sphere])
				
				# Rotate spheres
				#sphereSet[sphere].rotation_euler = (Euler((stepTime*simData[step][3 + properties*sphere],stepTime*simData[step][4 + properties*sphere],stepTime*simData[step][5 + properties*sphere])).to_matrix() @ sphereSet[sphere].rotation_euler.to_matrix()).to_euler()
				
				# Keyframe spheres
				sphereSet[sphere].keyframe_insert(data_path="location",index=-1)
				
				sphereSet[sphere].keyframe_insert(data_path="rotation_euler",index=-1)
			frameNum += 1

bpy.data.scenes[0].frame_end = frameNum
bpy.data.scenes[0].frame_start = 0 