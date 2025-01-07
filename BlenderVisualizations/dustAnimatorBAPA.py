from __future__ import division
import bpy
import numpy as np
from mathutils import *
from math import *
import fnmatch
import os

import sys
sys.path.append('/home/kolanzl/.local/lib/python3.11/site-packages')
import h5py


def get_filename(path,fileindex,relax = False):
    rel = ""
    if relax:
        rel = "RELAX"
        print("RELAX JOB")

        
    if os.path.exists(path):
        for file in os.listdir(path):
            if file.endswith(f"{rel}simData.csv"):
                filesplit = file.split("_")
                if len(filesplit) > 2: #Job Guidos way of naming
                    if fileindex == 0:
                        if not filesplit[1].isnumeric():
                            return file
                    else:
                        if filesplit[0] == str(fileindex) and filesplit[1].isnumeric():
                            return file
                else: #Lucas Kolanz way of naming
                    if filesplit[0] == str(fileindex):
                        return file
            elif file.startswith(str(fileindex)) and file.endswith(f"{rel}data.h5"):
                return file
    return -1
            
    

def get_simData_and_consts(path,fileindex,relax = False):
    numSpheres = -1
    steps = -1
    
    filename = get_filename(path,fileindex,relax)
    print(f"filename from get_filename: {filename}")
    
    if filename == -1:
        return [-1,-1,-1,-1]
    
    elif  filename.endswith("simData.csv"):
        
        try:
            print("fullpath: "+path+filename)
            simData = np.loadtxt(path+filename,dtype=float,delimiter=',',skiprows = 1)
        except Exception as e:
            print("ERROR CAUGHT in folder: {}".format(path+filename))
            print(e)
        #    simData = np.array([last_line.split(',')],dtype=np.float64)
#       print("fullpath")
#       print(path+file)
        print("DATA")
        print(simData)
        #simData = np.array([simData]) # Uncomment this line for single timestep data with no headers
        #simData.T # Uncomment this line for single timestep data with no headers
        steps = len(simData)
        print("steps: ",steps)
        constants = np.loadtxt(path + filename.replace("simData.csv", "constants.csv"),dtype=float,delimiter=',')
#       if fileindex == 0:
#           if filename.find("_") == filename.rfind("_"):
#           else:
#               constants = np.genfromtxt(path + filename[1:] + "constants.csv",dtype=float,delimiter=',')
#       else:
#           if "SpaceLab_branch" in path.split("/"):
#               constants = np.genfromtxt(path + str(fileindex) + filename[1:] + "constants.csv",dtype=float,delimiter=',')
#           else:
#               print("==================")
#               print(path + str(fileindex) + filename + "constants.csv")
#               print("==================")
#               constants = np.genfromtxt(path + str(fileindex) + filename + "constants.csv",dtype=float,delimiter=',')
        numSpheres = len(constants)
    elif filename.endswith(".h5"):
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
            if steps != 51:
                print("=========================================================")
                print(f"ERROR: there were {steps} writes in sim {filename}")
                print("=========================================================")
                    
    return [simData,constants,numSpheres,steps]
    

# Delete old stuff first:
foo_objs = [obj for obj in bpy.context.scene.objects if fnmatch.fnmatchcase(obj.name, "*phere*")]

for obj in foo_objs:
    bpy.data.objects.remove(obj, do_unlink = True)

foo_objs = [obj for obj in bpy.context.scene.objects if fnmatch.fnmatchcase(obj.name, "Mball*")]

for obj in foo_objs:
    bpy.data.objects.remove(obj, do_unlink = True)
    
foo_objs = [obj for obj in bpy.context.scene.objects if fnmatch.fnmatchcase(obj.name, "ProjectileArrow*")]

for obj in foo_objs:
    bpy.data.objects.remove(obj, do_unlink = True)
    

stepSkip = 1
stepTime = 1e-5
properties = 11
scaleUp = 1e5
frameNum = 0




path = '/media/kolanzl/easystore/SpaceLab_data/jobs/TESTBCCATWO7/N_256/T_3/'
path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA0/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA1/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA2/N_300/T_1000/' #good
path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA3/N_300/T_1000/' #bad cool ring
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA4/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA5/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA6/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA7/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA9/N_300/T_1000/' #bad
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA8/N_300/T_1000/' #good linear?
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPAloo0/N_300/T_1000/'


simStart = 30
simEnd = 300

#csv = False
#filename = ''  
just_last_line = False
job_group_suffix = path.split('/')[-4].split('_')[-1]
job_group_suffix = [i for i in job_group_suffix if not i.isnumeric()]
job_group_suffix = "".join(job_group_suffix)

rel = False
if job_group_suffix == "relax":
    rel = True

#filename = str(simStart)+"_simData.csv"

simData,constants,numSpheres,steps = get_simData_and_consts(path,simStart,relax=rel)

if isinstance(simData,int) and simData == -1:
    print(f"No data found for folder: {path}")
    raise ValueError("simData is -1 indicating an error in get_simData_and_consts")



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
    
#    sphereSet[sphere].hide_viewport = True
#    sphereSet[sphere].hide_render = True
#    sphereSet[sphere].hide_set(True)
    
    
#print(sphereSet)
#print("HERER")
#print(np.log2(simStart))
#print(np.log2(simEnd))
sims = [simStart*i for i in range(1,int(simEnd/simStart)+1)]
print(sims)
for i,sim in enumerate(sims):
    if frameNum > 1:
#        print(simData.shape)        
        simData,constants,numSpheres,steps = get_simData_and_consts(path,sim,rel)
#       print(sim)
#       print(simData.shape)
        print(f"total spheres in sim {i}: {numSpheres}")
        
        # Instanciaten the new particles and the end of file:
        for newball in range(sims[i-1],sims[i]):
#            print(newball)
            sphereSet.append(bpy.data.objects.new("Mball." + str(newball),sphereMesh))
#            print(len(sphereSet))
#            sphere += 1
    #       print(f"HERERE: {numSpheres-1}")
            bpy.context.scene.collection.objects.link(sphereSet[newball]) # link the object to the scene collection
            sphereSet[newball].scale = (scaleUp*constants[newball,0],scaleUp*constants[newball,0],scaleUp*constants[newball,0])
            sphereSet[newball].location = (scaleUp*simData[0][0 + properties*(newball)],scaleUp*simData[0][1 + properties*(newball)],scaleUp*simData[0][2 + properties*(newball)]) 
            sphereSet[newball].rotation_mode = "XYZ"
            
            sphereSet[newball].hide_viewport = True
            sphereSet[newball].hide_render = True
            sphereSet[newball].hide_set(True)
            sphereSet[newball].keyframe_insert(data_path="hide_viewport",index=-1)
            sphereSet[newball].keyframe_insert(data_path="hide_render",index=-1)
#            sphereSet[newball].keyframe_insert(data_path="hide_set",index=-1)

    
    if just_last_line:
        show_steps = [steps-1]
    else:
        show_steps = list(range(steps))
    for step in show_steps:
#        print("step: ",step) 
        if step % steps / 10 == 0:
            print(str(step/steps*100)+'%')
        if step % stepSkip == 0:
            bpy.context.scene.frame_set(frameNum)
            for sphere in range(numSpheres):
                
                if sphere < sim:
                    sphereSet[sphere].hide_viewport = False
                    sphereSet[sphere].hide_render = False
                    sphereSet[sphere].hide_set(False)
                    sphereSet[sphere].keyframe_insert(data_path="hide_viewport",index=-1)
                    sphereSet[sphere].keyframe_insert(data_path="hide_render",index=-1)
#                    sphereSet[sphere].keyframe_insert(data_path="hide_set",index=-1)
                else:
                    sphereSet[sphere].hide_viewport = True
                    sphereSet[sphere].hide_render = True
                    sphereSet[sphere].hide_set(True)
                    sphereSet[sphere].keyframe_insert(data_path="hide_viewport",index=-1)
                    sphereSet[sphere].keyframe_insert(data_path="hide_render",index=-1)
#                    sphereSet[sphere].keyframe_insert(data_path="hide_set",index=-1)
                
                # Move spheres
                sphereSet[sphere].location = (scaleUp*simData[step][0 + properties*sphere],scaleUp*simData[step][1 + properties*sphere],scaleUp*simData[step][2 + properties*sphere])
                
                # Rotate spheres
                #sphereSet[sphere].rotation_euler = (Euler((stepTime*simData[step][3 + properties*sphere],stepTime*simData[step][4 + properties*sphere],stepTime*simData[step][5 + properties*sphere])).to_matrix() @ sphereSet[sphere].rotation_euler.to_matrix()).to_euler()
                
                # Keyframe spheres
                sphereSet[sphere].keyframe_insert(data_path="location",index=-1)
                
                sphereSet[sphere].keyframe_insert(data_path="rotation_euler",index=-1)
                
#               if frameNum == 404
            frameNum += 1

bpy.data.scenes[0].frame_end = frameNum
bpy.data.scenes[0].frame_start = 0 