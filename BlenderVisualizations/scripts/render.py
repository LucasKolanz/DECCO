from __future__ import division

import bpy
import numpy as np
import os
import sys
import json
import fnmatch
import linecache
import mathutils
from math import *
import importlib

# -----------------------------------------------------------------------------
# Resolve project paths cleanly
# -----------------------------------------------------------------------------

# Folder containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

# Absolute path to SpaceLab project root
# (You can keep this if you're always running on the same machine)
project_path = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab/"

# Add utilities and visualizations folders
extra_paths = [
    os.path.join(project_path, "utilities"),
    os.path.join(project_path, "BlenderVisualizations"),
    "/home/kolanzl/.local/lib/python3.11/site-packages",   # optional
]

for p in extra_paths:
    if p not in sys.path:
        sys.path.append(p)

# -----------------------------------------------------------------------------
# Import SpaceLab modules (force clean reload for Blender)
# -----------------------------------------------------------------------------
import utils as u
import metricVisualization
importlib.reload(metricVisualization)
import metricVisualization as av



with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)

data_directory = input_json["data_directory"]
    

def get_filename(path,fileindex,relax = False):
    rel = ""
    if relax:
        rel = "RELAX"
        
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

        print("DATA")
        print(simData)

        steps = len(simData)
        print("steps: ",steps)
        constants = np.loadtxt(path + filename.replace("simData.csv", "constants.csv"),dtype=float,delimiter=',')

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
                print(f"WARNING: there were {steps} writes in sim {str(fileindex)+filename}")
                print("=========================================================")
        
                
    return [simData,constants,numSpheres,steps]
    

def place_point_light(relative_position):
    """
    Places a point light source at an XYZ position relative to the current camera.

    Args:
    - relative_position: A tuple representing the relative position of the point light
                         source with respect to the camera (X, Y, Z). The relative z axis
                         goes from the aggregate to the camera. The relative x points right 
                         of the camera and the relative y points up.
    """
    # Get the active camera
    camera = bpy.context.scene.camera
    if camera is None or camera.type != 'CAMERA':
        print("No active camera found!")
        return

    # Get the camera location
    camera_location = camera.location
    camera_rotation = camera.rotation_euler.to_matrix().to_4x4()

    # Calculate the relative position in world coordinates
    relative_position_world = camera_rotation @ mathutils.Vector(relative_position)

    # Calculate the point light source's position
    light_position = (
        camera_location[0] + relative_position_world.x,
        camera_location[1] + relative_position_world.y,
        camera_location[2] + relative_position_world.z
    )

    # Add a point light to the scene
    bpy.ops.object.light_add(type='POINT', location=light_position)
    point_light = bpy.context.object
    point_light.data.energy = 100000#15000#70000#100000  # Adjust light intensity as needed
    point_light.data.use_shadow = True

def render_scene(picName=""):
    
    set_background()


    if len(picName) == 0:
        picName = "testRender.png"

    scene = bpy.context.scene
    render = scene.render

    # ---------------------------------------------------------
    # Make sure filepath ends with .png BEFORE setting format
    # ---------------------------------------------------------
    render.use_file_extension = True

    if not picName.lower().endswith(".png"):
        picName += ".png"

    render.filepath = picName

    render.image_settings.file_format = 'PNG'
    render.image_settings.color_mode = 'RGB'
    render.image_settings.color_depth = '8'

    # ---------------------------------------------------------
    # Camera check
    # ---------------------------------------------------------
    cams = [obj for obj in bpy.data.objects if obj.type == 'CAMERA']
    if len(cams) != 1:
        delete_cameras()
        add_camera()

    # ---------------------------------------------------------
    # Render the frame
    # ---------------------------------------------------------
    bpy.ops.render.render(write_still=True)



def add_camera(camera_angle='default',shift=(0,0),focal_length=30.0):
    
    camera_dist = 60

     # Create a new camera
    bpy.ops.object.camera_add(location=(0, 0, 0))
    
    # Set the newly created camera as the active camera
    new_camera = bpy.context.object
    bpy.context.scene.camera = new_camera

    # Align the camera with the current view
    bpy.ops.view3d.camera_to_view_selected()

    # Set the camera's focal length to 1mm
    new_camera.data.lens = focal_length#30.0
    
    new_camera.data.shift_x = shift[0]
    new_camera.data.shift_y = shift[1]
    
    if camera_angle == 'default':
        # Move the camera back along the view vector
        #        mat = new_camera.matrix_world

        # Get the 3D view area
        area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')

        # Get the region data (this contains the view rotation)
        region_data = next(space for space in area.spaces if space.type == 'VIEW_3D').region_3d

        # Get the view direction
        view_direction = region_data.view_rotation @ mathutils.Vector((0, 0, -1))

        new_camera.rotation_euler = view_direction
        new_camera.location -= view_direction * camera_dist# Move the camera backwards by 100 meters
    else:

        camera_euler_angle = mathutils.Euler(camera_angle, 'XYZ')
        new_camera.rotation_euler = camera_angle
        location_unit_vec = camera_euler_angle.to_matrix() @ mathutils.Vector((0, 0, 1))
            # Move the camera backwards by 500 metersp
        new_camera.location = location_unit_vec*camera_dist


    place_point_light((13,0,-15))
    
    # Set the camera's display size
#    camera.data.display_size = 0.5
    new_camera.data.clip_end = 10000     
        
    print("new_camera.rotation_euler: {}".format(new_camera.rotation_euler))
    print("new_camera.location: {}".format(new_camera.location))
    return new_camera   

def perpendicular_axis(obj1, obj2):
    # Compute the vector connecting the two objects
    connecting_vector = obj2.location - obj1.location
    
    # Compute the midpoint of the two objects
    midpoint = (obj1.location + obj2.location) / 2
    
    # Find a vector that is not parallel to the connecting vector
    # A simple way is to use the world's up vector (0, 0, 1) unless the connecting vector is already vertical
    if connecting_vector.xy.length > 0.0001:  # Check if the vector is not vertical
        non_parallel_vector = mathutils.Vector((0, 0, 1))
    else:
        non_parallel_vector = mathutils.Vector((1, 0, 0))
    
    # Compute the cross product to get the perpendicular axis
    perpendicular = connecting_vector.cross(non_parallel_vector)
    perpendicular.normalize()  # Make it a unit vector
    
    return midpoint, perpendicular

def set_background():
    bpy.context.scene.view_settings.view_transform = 'Standard'

    scene = bpy.context.scene

    # Ensure a world exists
    if not scene.world:
        scene.world = bpy.data.worlds.new("World")

    world = scene.world
    world.use_nodes = True

    # Clear all nodes (important!)
    nt = world.node_tree
    nt.nodes.clear()

    # Create background and output nodes
    bg = nt.nodes.new("ShaderNodeBackground")
    out = nt.nodes.new("ShaderNodeOutputWorld")

    # Position nodes (optional, for cleanliness)
    bg.location = (0, 0)
    out.location = (200, 0)

    # Set white background
    bg.inputs[0].default_value = (1.0, 1.0, 1.0, 1.0)

    # Link them
    nt.links.new(bg.outputs["Background"], out.inputs["Surface"])

#    world = bpy.context.scene.world
#    world.use_nodes = True
#    bg = world.node_tree.nodes.get("Background")
#    bg.inputs[0].default_value = (1,1,1,1)

    
def add_sun():
    # Create a new sun light
    bpy.ops.object.light_add(type='SUN', align='WORLD', location=(0, 0, 0))

    # Optionally, you can set some properties for the sun light
    light = bpy.context.active_object
    light.rotation_euler = (0, radians(90), 0)

#    light.data.energy = 3  # Adjust the energy (brightness) of the sun
#    light.data.angle = 0.1  # Adjust the angle (size) of the sun, which affects the softness of shadows

#delete suns
def delete_cameras():
    for camera in bpy.data.cameras:
        bpy.data.cameras.remove(camera)
        
def delete_textures():
    for texture in bpy.data.textures:
        bpy.data.textures.remove(texture)
#delete suns
def delete_lights():
    for light in bpy.data.lights:
        bpy.data.lights.remove(light)
    
# Delete materials
def delete_materials():
    for material in bpy.data.materials:
        bpy.data.materials.remove(material)
        
def delete_collections():
    for collection in bpy.data.collections:
        if collection.name != "Master Collection":
            bpy.data.collections.remove(collection)

# Delete images
def delete_images():
    images_to_delete = [img for img in bpy.data.images]# if img.name.startswith("Image_Material")]
    for img in images_to_delete:
        bpy.data.images.remove(img)
        
def delete_objects():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()
       
def delete_meshes():
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)
    # Create a list of all mesh objects in the scene
#    mesh_objects = [obj for obj in bpy.context.scene.objects if obj.type == 'MESH']

#    # Unlink and delete each mesh object
#    for obj in mesh_objects:
#        bpy.context.scene.objects.unlink(obj)
#        bpy.data.objects.remove(obj)

#    # Clean up any orphaned mesh data
#    for mesh in bpy.data.meshes:
#        if mesh.users == 0:
#            bpy.data.meshes.remove(mesh)
#        
def clear_shaders():
    for obj in bpy.context.scene.objects:
        if obj.type == 'MESH' and obj.data.materials:
            for material in obj.data.materials:
                if material.use_nodes:
                    # Clear shader nodes from the material
                    material.node_tree.nodes.clear()
          
def delete_linestyles():
    for linestyle in bpy.data.linestyles:
        bpy.data.linestyles.remove(linestyle)
        
    try:
        view_layer = bpy.context.scene.view_layers["View Layer"]  # Replace "View Layer" with your view layer's name if different
    except KeyError:
        print("No linestyles to delete.")
        return
    # Iterate over the line sets and remove them
    while view_layer.freestyle_settings.linesets:
        view_layer.freestyle_settings.linesets.remove(view_layer.freestyle_settings.linesets[0])
        
              
def delete_actions():
    for act in bpy.data.actions:
        bpy.data.actions.remove(act)
        
        
        
# Create a new material with shader nodes
def create_black_material(num,opacity=1.0):
    material_name = "Black_Material{}".format(num)

    # Create a new material
    material = bpy.data.materials.new(name=material_name)

    # Use nodes for the material
    material.use_nodes = True
    shader_tree = material.node_tree

    for node in shader_tree.nodes:
        shader_tree.nodes.remove(node)

    # Add Principled BSDF shader node
    principled_bsdf = shader_tree.nodes.new(type='ShaderNodeBsdfPrincipled')
    principled_bsdf.location = (0, 0)
    
    # Set the base color to black
    howblack = 0.05
    principled_bsdf.inputs['Base Color'].default_value = (howblack,howblack,howblack,opacity)  # RGBA values, all set to 0 for black

    # Add Material Output node
    material_output = shader_tree.nodes.new(type='ShaderNodeOutputMaterial')
    material_output.location = (300, 0)
    shader_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])
    
    material.blend_method = 'HASHED'
#    material.shadow_method = 'HASHED'
    material.use_screen_refraction = True

    return material

# Create a new material with shader nodes
def create_material(red,green,blue):
    material_name = "Material"

    # Create a new material
    material = bpy.data.materials.new(name=material_name)

    # Use nodes for the material
    material.use_nodes = True
    shader_tree = material.node_tree

    for node in shader_tree.nodes:
        shader_tree.nodes.remove(node)

    # Add Principled BSDF shader node
    principled_bsdf = shader_tree.nodes.new(type='ShaderNodeBsdfPrincipled')
    principled_bsdf.location = (0, 0)
    
    # Set the base color to black
#    howblack = 0.05
    principled_bsdf.inputs['Base Color'].default_value = (red,green,blue,1)  # RGBA values, all set to 0 for black

    # Add Material Output node
    material_output = shader_tree.nodes.new(type='ShaderNodeOutputMaterial')
    material_output.location = (300, 0)
    shader_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])
    
    material.blend_method = 'HASHED'
#    material.shadow_method = 'HASHED'
    material.use_screen_refraction = True

    return material

def NthOpacity(N,opacity):
    index = np.argpartition(opacity,-N)[-N]
    return opacity[index]

def clearScene():
    # -------------------------------------------------------------------------
    # Clear scene
    # -------------------------------------------------------------------------
    print('deleting objects')
    delete_objects()

    print('delete lights')
    delete_lights()

    print('clear shaders')
    clear_shaders()

    print('deleting materials')
    delete_materials()

    print('delete actions')
    delete_actions()

    print('delete meshes')
    delete_meshes()

    print('deleting textures')
    delete_textures()

    print('deleting linestyles')
    delete_linestyles()

def makeFrames(path, camera_angle="default", shift=(0, 0), focal_length=30.0,
               addCameraAndLighting=False, addBalls=False, render=False, visualFile=None, visual_name = ""):

    # -------------------------------------------------------------------------
    # Basic rendering setup
    # -------------------------------------------------------------------------
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.transparent_max_bounces = 50
    bpy.context.scene.cycles.adaptive_threshold = 0.075

   

    # -------------------------------------------------------------------------
    # Camera & lighting
    # -------------------------------------------------------------------------
    if addCameraAndLighting:
        set_background()
        
        add_sun()

        print('deleting cameras')
        delete_cameras()

        camera = add_camera(camera_angle, shift, focal_length)

    # -------------------------------------------------------------------------
    # Simulation metadata
    # -------------------------------------------------------------------------
    stepSkip = 1
    stepTime = 1e-5
    properties = 11
    scaleUp = 1e5
    frameNum = 0

    A_index = -4
    N_index = -3
    T_index = -2

    print(path)

    # Extract group
    group = path.split("/")[A_index]
    print(f'group: {group}')

    job_group = group.split('_')[0]
    attempt = int(group.split("_")[-1])

    # Extract N
    N = int(path.split("/")[N_index].split("_")[-1])
    print(f"N: {N}")

    # Extract T
    Temp = int(path.split("/")[T_index].split("_")[-1])
    print(f"Temp: {Temp}")

    # Determine final timestep index
    simStart = u.find_max_index(path)

    # -------------------------------------------------------------------------
    # Load simulation data
    # -------------------------------------------------------------------------
    if addBalls:
        print(path)
        print(simStart)
        simData, constants, numSpheres, steps = get_simData_and_consts(
            path, simStart, relax=False
        )
        print(np.array(simData).shape)

        if isinstance(simData, int) and simData == -1:
            print(f"No data found for folder: {path}")
            raise ValueError("simData is -1 indicating an error in get_simData_and_consts")

        try:
            simData.shape[1]
            simData = simData[-1]
        except Exception:
            pass

        print("simdata shape: " + str(simData.shape))
        print("constData shape: " + str(constants.shape))

        # -------------------------------------------------------------------------
        # Create base sphere mesh
        # -------------------------------------------------------------------------
        print('create mesh')
        bpy.ops.mesh.primitive_uv_sphere_add(location=(0, 0, 0), radius=1)
        obj = bpy.context.object
        sphereMesh = obj.data
        bpy.data.objects.remove(obj)

        print("num spheres: ", numSpheres)

        # -------------------------------------------------------------------------
        # Create collection for spheres
        # -------------------------------------------------------------------------
        sphere_collection_name = "Spheres"

        if sphere_collection_name in bpy.data.collections:
            sphere_collection = bpy.data.collections[sphere_collection_name]
        else:
            sphere_collection = bpy.data.collections.new(sphere_collection_name)
            bpy.context.scene.collection.children.link(sphere_collection)

        # -------------------------------------------------------------------------
        # Material
        # -------------------------------------------------------------------------
        howBlack = 0.05
        material = create_material(howBlack, howBlack, howBlack)
        material.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = 1.0

        sphereSet = []

        print("Instantiate spheres")

        # -------------------------------------------------------------------------
        # Instantiate sphere objects
        # -------------------------------------------------------------------------
        if len(sphereSet) == 0:
            for sphere_num in range(numSpheres):
                ballname = f"Mball.{sphere_num}"
                sphere = bpy.data.objects.new(ballname, sphereMesh)

                # Link to sphere collection
                sphere_collection.objects.link(sphere)

                bpy.context.view_layer.objects.active = sphere

                if not sphere.material_slots:
                    bpy.ops.object.material_slot_add()

                sphere.material_slots[0].link = 'OBJECT'
                sphereSet.append(sphere)

        print(simData.shape)

        # -------------------------------------------------------------------------
        # Apply scale, location, rotation to each sphere
        # -------------------------------------------------------------------------
        for sphere_num in range(numSpheres):
            sphere = sphereSet[sphere_num]
            sphere.material_slots[0].material = material

            sphere.scale = (scaleUp * constants[sphere_num, 0],) * 3

            sphere.location = (
                scaleUp * simData[properties * sphere_num + 0],
                scaleUp * simData[properties * sphere_num + 1],
                scaleUp * simData[properties * sphere_num + 2]
            )

            sphere.rotation_mode = "XYZ"

    # -------------------------------------------------------------------------
    # Optional rendering
    # -------------------------------------------------------------------------
    if render:
        render_scene(
            f"{data_directory}data/figures/aggRenders/metricVisuals/"
            f"visual_{visual_name}-{job_group}_a-{attempt}_N-{N}_T-{Temp}.png"
        )

    print('================================end================================')

    
    
    

paths = []
camera_dic = {}



######################################################################################################################
#Visualize all the metrics
paths.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_12/N_300/T_3/')
camera_dic[paths[0]] = (-1.5708, 0.0000, -2.0176)

#clearScene()
#av.addEquivalentEllipsoid(paths[-1])
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True,visual_name="Pabc")

#clearScene()
#av.addGyrationRadiusSphere(paths[-1])
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True,visual_name="PKBM")

#clearScene()
#av.addEnclosingEllipsoid(paths[-1])
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True,visual_name="Pfee")

#clearScene()
av.addEnclosingSphere(paths[-1])
makeFrames(paths[-1],camera_dic[paths[-1]],render=True,addBalls=False,addCameraAndLighting=False,visual_name="Pfes")

#clearScene()
#av.addConvexHull(paths[-1])
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True,visual_name="Pch")

#clearScene()
#av.addGCSSphere(paths[-1])
#av.addShadowGrid(paths[-1],grid_offset=10)
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True,addCameraAndLighting=False,addBalls=False,visual_name="gcs")

######################################################################################################################
