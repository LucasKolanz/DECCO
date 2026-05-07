from __future__ import division
import bpy
import numpy as np
from mathutils import *
from math import *
import fnmatch

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

#path = "C:/Users/milin/Desktop/VSOUT/1000-0/150_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
#path = "C:/Users/milin/Desktop/VSOUT/3-0/220_2_rho2.25_k5e-03_Ha5e-12_dt5e-10_"
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

path = '/home/lucas/Desktop/Research/SpaceLabTesting/SpaceLab/output/dust_collision/'
#path = '/home/lucas/Desktop/Research/SpaceLabTesting/SpaceLab/ColliderSingleCore/'
#filename = '_20_R1e-04_v4e-01_cor0.63_mu0.1_rho2.25_k1e+01_Ha5e-12_dt3e-10_'
#filename = '1_19_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_'
filename = '19_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_19_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_T0.00_vBig4e+00_vSmall-4e+00_IP0._rho2.25_k4e+02_Ha5e-12_dt5e-11_'

simData = np.loadtxt(path + filename + "simData.csv",dtype=float,delimiter=',',skiprows = 1)
#simData = np.array([simData]) # Uncomment this line for single timestep data with no headers
#simData.T # Uncomment this line for single timestep data with no headers
steps = len(simData)
#print(steps)
constants = np.genfromtxt(path + filename + "constants.csv",dtype=float,delimiter=',')
numSpheres = len(constants)

sphereSet = []
actionSet = []

# Initial sphere mesh to be instanced:
bpy.ops.mesh.primitive_uv_sphere_add(location = (0,0,0), radius = 1)
bpy.data.objects[2].hide_set(True)
#bpy.ops.object.metaball_add(type='BALL', radius=2, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
obj = bpy.context.object # the currently selected object
#obj.data.resolution = .1
#obj.data.threshold = 1.7
sphereMesh = obj.data # retrieve the mesh
#bpy.data.objects.remove(obj) # remove the object



# Instanciate spheres:
for sphere in range(numSpheres):
    sphereSet.append(bpy.data.objects.new("Mball." + str(sphere),sphereMesh))
    bpy.context.scene.collection.objects.link(sphereSet[sphere]) # link the object to the scene collection
    sphereSet[sphere].scale = (scaleUp*constants[sphere,0],scaleUp*constants[sphere,0],scaleUp*constants[sphere,0])
    sphereSet[sphere].location = (scaleUp*simData[0][0 + properties*sphere],scaleUp*simData[0][1 + properties*sphere],scaleUp*simData[0][2 + properties*sphere])
    sphereSet[sphere].rotation_mode = "XYZ"

frameNum = 0
for step in range(steps):
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