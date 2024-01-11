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
     
#path = '/home/kolanzl/Desktop/SpaceLab/jobs/SETTEST/T_3/'
#path = '/home/lucas/Desktop/Research/SpaceLabTesting/SpaceLab/ColliderSingleCore/'
#filename = '_20_R1e-04_v4e-01_cor0.63_mu0.1_rho2.25_k1e+01_Ha5e-12_dt3e-10_'
#filename = '_2_R5e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_'
#sim = 1
#simData = np.loadtxt(path + str(sim) + filename + "simData.csv",dtype=float,delimiter=',',skiprows = 1)
#simData = np.array([simData]) # Uncomment this line for single timestep data with no headers
#simData.T # Uncomment this line for single timestep data with no headers
#steps = len(simData)
#print("steps: ",steps)
#constants = np.genfromtxt(path + str(sim) + filename + "constants.csv",dtype=float,delimiter=',')
#numSpheres = len(constants)
#print("num spheres: ",numSpheres)

sphereSet = []
actionSet = []

radius = 2

# Initial sphere mesh to be instanced:
#bpy.ops.mesh.primitive_uv_sphere_add(location = (1,1,1), radius = 1)
bpy.ops.object.metaball_add(type='BALL', radius=radius, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
#bpy.data.objects[2].hide_set(True)
obj = bpy.context.object # the currently selected object
#obj.data.resolution = .1
#obj.data.threshold = 1.7
sphereMesh = obj.data # retrieve the mesh
#bpy.data.objects.remove(obj) # remove the object

#xcm = 0
#ycm = 0
#zcm = 0
#mtot = 0
#for sphere in range(2,12):
#    xcm += simData[0][0 + properties*sphere] * constants[sphere,1]
#    ycm += simData[0][1 + properties*sphere] * constants[sphere,1]
#    zcm += simData[0][2 + properties*sphere] * constants[sphere,1]
#    mtot += constants[sphere,1]
#xcm /= mtot
#ycm /= mtot
#zcm /= mtot
    
# Instanciate spheres:
sphere = 0
sphereSet.append(bpy.data.objects.new("Mball." + str(sphere),sphereMesh))
bpy.context.scene.collection.objects.link(sphereSet[sphere]) # link the object to the scene collection
#sphereSet[sphere].scale = (scaleUp*constants[sphere,0],scaleUp*constants[sphere,0],scaleUp*constants[sphere,0])
sphereSet[sphere].location = (0,0,0)      
sphereSet[sphere].rotation_mode = "XYZ"

sphere = 1
sphereSet.append(bpy.data.objects.new("Mball." + str(sphere),sphereMesh))
bpy.context.scene.collection.objects.link(sphereSet[sphere]) # link the object to the scene collection
#sphereSet[sphere].scale = (scaleUp*constants[sphere,0],scaleUp*constants[sphere,0],scaleUp*constants[sphere,0])
sphereSet[sphere].location = (0,6,0)      
sphereSet[sphere].rotation_mode = "XYZ"
