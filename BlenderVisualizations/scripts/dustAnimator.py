from __future__ import division
import bpy
import numpy as np
from collections import deque
from mathutils import *
from math import *
import fnmatch
import os
import re

import sys
sys.path.append('/home/kolanzl/.local/lib/python3.11/site-packages')
#import h5py

# relative_path = "../../"
# #print("\n\n")
# #print(relative_path)
# relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
# #print(relative_path)
# project_path = os.path.abspath(relative_path) + '/'
# sys.path.append(project_path+"utilities/")
# ## sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
# import utils as u


def unwrap_periodic_positions(positions, box_size, box_min=None, recenter=True):
    """
    Unwrap an aggregate split across periodic boundaries.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Wrapped particle positions.
    box_size : float or array_like, shape (3,)
        Periodic box length. A scalar assumes a cubic box.
    box_min : float or array_like, shape (3,), optional
        Lower edge of the periodic box. Defaults to -box_size / 2.
        Use 0 if your box runs from 0 to box_size.
    recenter : bool, default=True
        If True, move the aggregate's center of mass to the box center.

    Returns
    -------
    unwrapped : ndarray, shape (N, 3)
        Contiguous particle positions. These may temporarily lie outside
        the original periodic box when recenter=False.
    """
    positions = np.asarray(positions, dtype=float)

    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")

    box_size = np.broadcast_to(
        np.asarray(box_size, dtype=float), 3
    ).copy()

    if np.any(box_size <= 0):
        raise ValueError("box_size must be positive")

    if box_min is None:
        box_min = -box_size / 2
    else:
        box_min = np.broadcast_to(
            np.asarray(box_min, dtype=float), 3
        ).copy()

    # First map everything into [0, L).
    wrapped = (positions - box_min) % box_size
    unwrapped = wrapped.copy()

    for axis in range(3):
        L = box_size[axis]
        order = np.argsort(wrapped[:, axis])
        sorted_x = wrapped[order, axis]

        # Gaps between adjacent particles, including the periodic gap.
        gaps = np.diff(np.concatenate((sorted_x, [sorted_x[0] + L])))
        largest_gap_index = np.argmax(gaps)

        # Put the periodic cut immediately after the largest empty gap.
        cut = sorted_x[(largest_gap_index + 1) % len(sorted_x)]

        unwrapped[:, axis] = (wrapped[:, axis] - cut) % L

    if recenter:
        aggregate_center = np.mean(unwrapped, axis=0)
        box_center = box_min + box_size / 2
        unwrapped += box_center - aggregate_center
    else:
        unwrapped += box_min

    return unwrapped

#finds the integer corresponding to a particular part of the directory structure
#could be T for temperature, N for number of particles, etc.
def value_from_directory(pattern,directory):
    if pattern != "job_group" and pattern != "attempt" and pattern != "a":
        # Matches integers (e.g., 42) or decimals (e.g., 42.5)
        val = re.search(rf"{re.escape(pattern)}_([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", directory)
        
        if val:
            num_str = val.group(1)
            # Keeps distinct types: float if a dot is present, otherwise int
            numVal = float(num_str) if '.' in num_str else int(num_str)
        else:
            numVal = None
            
        return numVal

    else:
        split_dir = directory.split("/")[::-1]
        split_again = [folder.split('_') for folder in split_dir]
        # print(split_again)

        group_index = -1
        for i in range(len(split_dir)):
            # if split_again[i][0] == "N":
            # print(split_again[i][0])
            # print(i)
            # print(split_again[i+1][0])
            # print(i+1)
            if split_again[i][0] == "N" and split_again[i+1][0] == "M":
                group_index = i+2
            elif split_again[i][0] == "N" and split_again[i+1][0] != "M":
                group_index = i+1
        if group_index < 0:
            print(f"Group index not found for {directory}.\nNow exiting...")
            exit(-1)

        if pattern == "job_group":
            return split_again[group_index][0]
        elif pattern == "attempt" or pattern == 'a':
            return int(split_again[group_index][1])
        else:
            print("This should never happen :))")
            exit(-1)

#returns all indices in a directory if {index}_checkpoint.txt exists, or 
#if the data file exits and timing.txt exists.
def get_all_indices(directory):
    indices = []
    is_finished = os.path.exists(directory+"timing.txt")
    print(is_finished)
    for file in os.listdir(directory):
        if file.split("_")[0].isnumeric():
            #if this index has checkpointed, or this is the data file and the simulation has finished
            if file.endswith("checkpoint.txt") or ((file.endswith("simData.csv") or (file.endswith("data.h5"))) and is_finished):
                index = int(file.split("_")[0])
                #make sure this index isn't already in indices
                if index not in indices:
                    indices.append(index)
    print(indices)
    return list(sorted(set(indices)))


def get_filename(path,fileindex,relax = False):
    rel = ""
    if relax:
        rel = "RELAX"
        print("RELAX JOB")

    print(f"path in get_filename:{path}")
        
    if os.path.exists(path):
        for file in os.listdir(path):
            if file.endswith(f"{rel}simData.csv"):
#                print(file)
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
    else:
        print(f"ERROR: Path does not exist: {path}")
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
            print(simData.shape)
#            if len(simData.shape):
#                sim
        except Exception as e:
            print("ERROR CAUGHT in folder: {}".format(path+filename))
            print(e)
        #    simData = np.array([last_line.split(',')],dtype=np.float64)
#		print("fullpath")
#		print(path+file)
#        print("DATA")
#        print(simData)
        #simData = np.array([simData]) # Uncomment this line for single timestep data with no headers
        #simData.T # Uncomment this line for single timestep data with no headers
        steps = len(simData)
        print("steps: ",steps)
        constants = np.loadtxt(path + filename.replace("simData.csv", "constants.csv"),dtype=float,delimiter=',')
#		if fileindex == 0:
#			if filename.find("_") == filename.rfind("_"):
#			else:
#				constants = np.genfromtxt(path + filename[1:] + "constants.csv",dtype=float,delimiter=',')
#		else:
#			if "SpaceLab_branch" in path.split("/"):
#				constants = np.genfromtxt(path + str(fileindex) + filename[1:] + "constants.csv",dtype=float,delimiter=',')
#			else:
#				print("==================")
#				print(path + str(fileindex) + filename + "constants.csv")
#				print("==================")
#				constants = np.genfromtxt(path + str(fileindex) + filename + "constants.csv",dtype=float,delimiter=',')
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
       
    # print("HERE")
    # print(simData)             
    return [simData,constants,numSpheres,steps]



def pos_m_r(simData,constants):
    pos = simData.reshape(simData.shape[0], -1, 11)[:, :, :3]#[step,ball,xyz]
    m = constants[:,1]
    r = constants[:,0]

    return pos,m,r


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
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA3/N_300/T_1000/' #bad cool ring
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA4/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA5/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA6/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA7/N_300/T_1000/' #good
#path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA9/N_300/T_1000/' #bad
path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA8/N_300/T_1000/' #good linear?
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/SeqStickLognorm_1/N_300/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/SeqStickLognormrelax_1/N_300/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/test/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/testest_0/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/constrollingfric0/N_300/T_300/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsNovus/constrelax_4/N_300/T_30/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/JKRTestest_0/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/spinTest0/N_3/T_3/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/WELDTEST_1/M_10/N_30/T_1000/'
path = '/home/kolanzl/Desktop/Visualize/V3/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_branch/SpaceLab_data/jobs/bigbox_5/N_1000/B_3.000e-03/T_3/'
path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/weldVerificationCollision_6/'
#9,12,

simStart = 300
simEnd = 300
just_last_line = False
bigbox = False

startline = 0
endline = -1
stepline = 1


if endline == -1:
    with open(path + f"{simStart}_simData.csv", "rb") as f:
        endline =  sum(chunk.count(b"\n") for chunk in iter(lambda: f.read(1024 * 1024), b"")) -1


#csv = False
#filename = ''  
job_group_suffix = path.split('/')[-4].split('_')[-1]
job_group_suffix = [i for i in job_group_suffix if not i.isnumeric()]
job_group_suffix = "".join(job_group_suffix)
 
rel = ("relax" in path)
if "relax" in path:
    rel = True


#filename = str(simStart)+"_simData.csv"

sims = get_all_indices(path)
#sims = [3]

sims = [i for i in sims if i >= simStart and i <= simEnd]
print(f"sims: {sims}")


simData,constants,numSpheres,steps = get_simData_and_consts(path,sims[0],relax=rel)
if len(simData.shape) == 1:
    simData = simData[np.newaxis,:]
if isinstance(simData,int) and simData == -1:
    print(f"No data found for folder: {path}")
    raise ValueError("simData is -1 indicating an error in get_simData_and_consts")
# print(simData)
# print(simData.shape)
pos,m,r = pos_m_r(simData,constants) 

if bigbox or "bigbox" in path:
    boxsize = value_from_directory("B",path)
    for i,p in enumerate(pos):
        pos[i] = unwrap_periodic_positions(p,boxsize)  


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
    sphereSet[sphere].scale = (scaleUp*r[sphere],scaleUp*r[sphere],scaleUp*r[sphere])
    sphereSet[sphere].location = (scaleUp*pos[0,sphere,0],scaleUp*pos[0,sphere,1],scaleUp*pos[0,sphere,2])
    # sphereSet[sphere].scale = (scaleUp*constants[sphere,0],scaleUp*constants[sphere,0],scaleUp*constants[sphere,0])
    # sphereSet[sphere].location = (scaleUp*simData[0][0 + properties*sphere],scaleUp*simData[0][1 + properties*sphere],scaleUp*simData[0][2 + properties*sphere])
#    print(simData[0][0 + properties*sphere],simData[0][1 + properties*sphere],simData[0][2 + properties*sphere])
    sphereSet[sphere].rotation_mode = "XYZ"
    
#    sphereSet[sphere].hide_viewport = True
#    sphereSet[sphere].hide_render = True
#    sphereSet[sphere].hide_set(True)
    
    
#print(sphereSet)
#print("HERER")
#print(np.log2(simStart))
#print(np.log2(simEnd))
#sims = [simStart*i for i in range(1,int(simEnd/simStart)+1)]
#if len(sims) = 1:
#    sims = 
for i,sim in enumerate(sims):
    print(f"getting sim index {sim}")
    if frameNum > 1:
#        print(simData.shape)        
        simData,constants,numSpheres,steps = get_simData_and_consts(path,sim,rel)
        pos,m,r = pos_m_r(simData,constants)
        if bigbox or "bigbox" in path:
            boxsize = value_from_directory("B",path)
            for i,p in enumerate(pos):
                pos[i] = unwrap_periodic_positions(p,boxsize)
        print(f"Steps: {steps}")
#	    print(simData.shape)
        print(f"total spheres in sim {i}: {numSpheres}")
        
        # Instanciaten the new particles and the end of file:
        for newball in range(sims[i-1],sims[i]):
#            print(newball)
            sphereSet.append(bpy.data.objects.new("Mball." + str(newball),sphereMesh))
#            print(len(sphereSet))
#            sphere += 1
    #		print(f"HERERE: {numSpheres-1}")
            bpy.context.scene.collection.objects.link(sphereSet[newball]) # link the object to the scene collection
            sphereSet[newball].scale = (scaleUp*r[newball],scaleUp*r[newball],scaleUp*r[newball])
            sphereSet[newball].location = (scaleUp*pos[0,newball,0],scaleUp*pos[0,newball,1],scaleUp*pos[0,newball,2])
            # sphereSet[newball].scale = (scaleUp*constants[newball,0],scaleUp*constants[newball,0],scaleUp*constants[newball,0])
            # sphereSet[newball].location = (scaleUp*simData[0][0 + properties*(newball)],scaleUp*simData[0][1 + properties*(newball)],scaleUp*simData[0][2 + properties*(newball)]) 
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
        # show_steps = list(range(steps))
        show_steps = list(range(startline,endline,stepline))
    for step in show_steps:
        if step%stepSkip == 0:
            print("step: ",step) 
            nan_data = (np.isnan(simData[step][0 + properties*0]) or np.isnan(simData[step][1 + properties*0]) or np.isnan(simData[step][2 + properties*0]))
            if not nan_data:
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
                        sphereSet[sphere].location = (scaleUp*pos[step,sphere,0],scaleUp*pos[step,sphere,1],scaleUp*pos[step,sphere,2])
                        # sphereSet[sphere].location = (scaleUp*simData[step][0 + properties*sphere],scaleUp*simData[step][1 + properties*sphere],scaleUp*simData[step][2 + properties*sphere])
                        
                        # Rotate spheres
                        # sphereSet[sphere].rotation_euler = (Euler((stepTime*simData[step][3 + properties*sphere],stepTime*simData[step][4 + properties*sphere],stepTime*simData[step][5 + properties*sphere])).to_matrix() @ sphereSet[sphere].rotation_euler.to_matrix()).to_euler()
                        
                        # Keyframe spheres
                        sphereSet[sphere].keyframe_insert(data_path="location",index=-1)
                        
                        sphereSet[sphere].keyframe_insert(data_path="rotation_euler",index=-1)
                        
    #				if frameNum == 404
                frameNum += 1

bpy.data.scenes[0].frame_end = frameNum
bpy.data.scenes[0].frame_start = 0  