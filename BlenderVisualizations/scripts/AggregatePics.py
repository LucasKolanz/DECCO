from __future__ import division

import os
import sys

project_path = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab/"

extra_paths = [
    "/home/kolanzl/.local/lib/python3.13/site-packages",
    os.path.join(project_path, "utilities"),
    os.path.join(project_path, "BlenderVisualizations"),
]

for p in reversed(extra_paths):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)

sys.path = list(dict.fromkeys(sys.path))

import bpy
import numpy as np
import h5py
import mathutils
from mathutils import *
from math import *
import fnmatch
import json
import linecache

import utils as u


with open(project_path+"default_files/default_input.json",'r') as fp:
    input_json = json.load(fp)

data_directory = input_json["data_directory"]
    

def hex_to_rgb(hex_color):
    """
    Convert a hex color code to a normalezed RGB tuple.

    """

    hex_color = hex_color.strip().lstrip("#")

    if len(hex_color) != 6:
        raise ValueError("Hex color must have exactly 6 hexadecimal characters")

    nonnormalized = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return tuple(c*1.0/255.0 for c in nonnormalized)

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
    point_light.data.energy = 175000#100000#15000#70000#100000  # Adjust light intensity as needed
    point_light.data.use_shadow = True

def render_scene(picName=""):

    if len(picName) == 0:
        picName = "testRender.png"

    # Ensure there's only one camera in the scene
    if len([obj for obj in bpy.data.objects if obj.type == 'CAMERA']) > 1:
        delete_cameras()
        add_camera()
    elif len([obj for obj in bpy.data.objects if obj.type == 'CAMERA']) == 0:
        add_camera()
    #    else:
    #        print("There's either no camera or more than one camera in the scene.")
    
    # Set render settings
    scene = bpy.context.scene
    render = scene.render
    render.image_settings.file_format = 'PNG'
    render.image_settings.color_mode = 'RGB'  # RGB if no transparency needed
    render.image_settings.color_depth = '8'   # '8' or '16'
    render.filepath = picName

    # Resolution (e.g., 4K)
#    render.resolution_x = 3840
#    render.resolution_y = 2160
#    render.resolution_percentage = 100

    # Use Cycles for high quality
#    scene.render.engine = 'CYCLES'
#    scene.cycles.device = 'GPU'  # or 'CPU'

    # Sampling and denoising
#    scene.cycles.samples = 512   # Increase for better quality
#    scene.cycles.use_denoising = False

    # Optional: transparent background
#    scene.render.film_transparent = True

    # Render and write to file
    bpy.ops.render.render(write_still=True)

def add_camera(camera_angle='default',shift=(0,0),focal_length=30.0,camera_dist=60):

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

        camera_euler_angle = Euler(camera_angle, 'XYZ')
        new_camera.rotation_euler = camera_angle
        location_unit_vec = camera_euler_angle.to_matrix() @ mathutils.Vector((0, 0, 1))
            # Move the camera backwards by 500 metersp
        new_camera.location = location_unit_vec*camera_dist


    # place_point_light((13,0,-15))
    
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
    # Ensure there's a world in the scene
    if not bpy.context.scene.world:
        bpy.context.scene.world = bpy.data.worlds.new("World")

    world = bpy.context.scene.world

    # Use nodes
    world.use_nodes = True

    # Get the background node
    bg_node = world.node_tree.nodes.get('Background')
    if not bg_node:
        bg_node = world.node_tree.nodes.new(type='ShaderNodeBackground')

    # Set the background color to white
    bg_node.inputs[0].default_value = (1, 1, 1, 1)  # RGBA for white
    
def add_sun(
    energy=1.0,
    camera=None,
    direction_camera_frame=(1, 0, -1)
):
    """
    Add a sun whose direction is defined relative to the camera frame.

    direction_camera_frame is the direction the light travels in the camera's
    local coordinates.

    Camera frame convention in Blender:
        local -Z = camera looks forward
        local +Y = camera up
        local +X = camera right
    """

    if camera is None:
        camera = bpy.context.scene.camera

    if camera is None:
        raise ValueError("No camera provided and scene has no active camera.")

    # Add sun
    bpy.ops.object.light_add(type='SUN', align='WORLD', location=(0, 0, 0))
    light = bpy.context.active_object

    # Set brightness
    light.data.energy = energy

    # Direction of incoming sunlight in camera-local coordinates
    local_dir = Vector(direction_camera_frame).normalized()

    # Convert that direction from camera frame to world frame
    world_dir = camera.matrix_world.to_3x3() @ local_dir
    world_dir.normalize()

    # A Blender light points along its local -Z axis.
    # Rotate the sun so its -Z axis points in world_dir.
    light.rotation_euler = world_dir.to_track_quat('-Z', 'Y').to_euler()

    return light


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


def makeFrames(path,camera_angle="default",shift=(0,0),focal_length=30.0,camera_dist=60,color=(),sun_energy=1.0,render=False):    
        
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.transparent_max_bounces = 50
    bpy.context.scene.cycles.adaptive_threshold = 0.075

#    bpy.context.space_data.clip_end = 50000

    print('deleting objects')
    delete_objects()    
    print('delete lights')
    delete_lights()
    print('clear shaders')
    clear_shaders()
    print('deleting materials')
    delete_materials()
#    print("delete images")
#    delete_images()
    print('delete actions')
    delete_actions()
    print('delete meshes')
    delete_meshes()
    print('deleting textures')
    delete_textures()
    print('deleting linestyles')
    delete_linestyles()

    
    print('deleting cameras')
    delete_cameras()
#        print("Camera info from doEverything: {}, {}, {}".format(camera_angle,camera_position,camera_name))
    camera=add_camera(camera_angle,shift,focal_length,camera_dist)
    set_background()
    add_sun(sun_energy,camera=camera) #This has to come after the camera is added because it is placed based on the camera

    stepSkip = 1
    stepTime = 1e-5
    properties = 11
    scaleUp = 1e5
    frameNum = 0
    
    A_index = -5
    # A_index = -4
#    M_index = -4
    N_index = -3
    T_index = -2
    
    print(path)

    # att = ""
    # group = path.split("/")[A_index]
    # print(f'group: {group}')
    # job_group = group.split('_')[0]
    # attempt = int(group.split("_")[-1])
    attempt = u.value_from_directory("attempt",path)
    job_group = u.value_from_directory("job_group",path)
    print(f"path: {path}")
    print(f"attempt: {attempt}")
    print(f"job_group: {job_group}")

#    M = int(path.split("/")[M_index].split("_")[-1])
#    print(f"M: {M}")
    M = u.value_from_directory("M",path)
    if M is None:
        M = 1
    print(f"M: {M}")


    # N = int(path.split("/")[N_index].split("_")[-1])
    N = u.value_from_directory("N",path)
    print(f"N: {N}")

    # Temp = int(path.split("/")[T_index].split("_")[-1])
    Temp = u.value_from_directory("T",path)
    print(f"Temp: {Temp}")

    
#    simStart = N-3
    simStart = u.find_max_index(path)
     
    
#    job_group = job_group.strip(str(attempt))
    
    #add the ellipsoid
#    addEllipse(path,N)
#    addCircle(path,N)
    
    print(path)
    print(simStart)
    simData,constants,numSpheres,steps = get_simData_and_consts(path,simStart,relax=False)
    print(np.array(simData).shape)

    if isinstance(simData,int) and simData == -1:
        print(f"No data found for folder: {path}")
        raise ValueError("simData is -1 indicating an error in get_simData_and_consts")

    try:
        simData.shape[1]
        simData = simData[-1]
    except IndexError:
        pass
    print("simdata shape: " + str(simData.shape))
    print("constData shape: " + str(constants.shape))


    ## Initial sphere mesh to be instanced:
    print('create mesh')
    bpy.ops.mesh.primitive_uv_sphere_add(
        location=(0, 0, 0),
        radius=1,
        # segments=8,     # longitude divisions
        # ring_count=4    # latitude divisions
        segments=128,     # longitude divisions
        ring_count=64    # latitude divisions
    )
    obj = bpy.context.object # the currently selected object
    sphereMesh = obj.data # retrieve the mesh
    bpy.data.objects.remove(obj) # remove the object



    print("num spheres: ",numSpheres)
    


    #make material
    material = create_material(color[0],color[1],color[2])
    #material = create_black_material(1)
    # temps = [3,10,30,100,300,1000]
    # temp = Temp
    # index = temps.index(temp)
    # blue_percent = 1.0-index*1.0/((len(temps)-1)*1.0) 
    # red_percent = index*1.0/((len(temps)-1)*1.0)
    # green_percent = 0
#    print(f"index: {index}")
#    print(f"blue_percent: {blue_percent}")
#    print(f"red_percent: {red_percent}")
    # material = create_material(red_percent,green_percent,blue_percent)
#    material = create_material(0.9686274509803922,0.4980392156862745,0) #OSU orange
#    howBlack = 0.05
#    material = create_material(howBlack,howBlack,howBlack)
#    material = create_material(0.075,0.075,1)
#    material = create_material(0,0.4,0) #green
#    offset = -0.15
#    material = create_material(0.4235+offset+0.3,0.2314+offset,0.6667+offset)
    material.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = 1.0
    sphereSet = []

    print('Instanciate spheres')
    # Instanciate spheres:
#        bpy.context.scene.frame_set(frameNum)
    if len(sphereSet) == 0:
        for sphere_num in range(numSpheres):
            ballname = f"Mball.{sphere_num}"
            sphere = bpy.data.objects.new(ballname,sphereMesh)
            bpy.context.scene.collection.objects.link(sphere) # link the object to the scene collection
            bpy.context.view_layer.objects.active = sphere
            if not sphere.material_slots:

                bpy.ops.object.material_slot_add()
#                bpy.ops.object.material_slot_add({'object': sphere})
                
            sphere.material_slots[0].link = 'OBJECT'
            sphereSet.append(sphere)

    # print(simData.shape)
    for sphere_num in range(numSpheres):
        # print(sphere_num)
        sphereSet[sphere_num].material_slots[0].material = material
#        materials[sphere_num].node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = opacityData[sphere_num]
        sphereSet[sphere_num].scale = (scaleUp*constants[sphere_num,0],scaleUp*constants[sphere_num,0],scaleUp*constants[sphere_num,0])
        sphereSet[sphere_num].location = (scaleUp*simData[properties*sphere_num+0],scaleUp*simData[properties*sphere_num+1],scaleUp*simData[properties*sphere_num+2])
        sphereSet[sphere_num].rotation_mode = "XYZ"
       
    # print(f"{data_directory}data/figures/aggRenders/BAPA/Coloredagg-{job_group}_a-{attempt}_N-{N}_T-{Temp}.png")
    if render: 
#        render_scene(f"{data_directory}data/figures/aggRenders/agg-{group}_a-{attempt}_N-{N}_seqStick.png")
#        render_scene(f"{data_directory}data/figures/aggRenders/agg-{group}_a-{attempt}_M-{M}_N-{N}_T-{Temp}_niether.png")
        # render_scene(f"{data_directory}data/figures/aggRenders/Coloredagg-{job_group}_a-{attempt}_N-{N}_T-{Temp}.png")
       # render_scene(f"{data_directory}data/figures/aggRenders/ColoredFragg-{job_group}_a-{attempt}_M-{M}_T-{Temp}.png")
       render_scene(f"{data_directory}data/figures/aggRenders/paper2/ColoredFragg-{job_group}_a-{attempt}_M-{M}_N-{N}_T-{Temp}.png")
       # render_scene(f"{data_directory}data/figures/aggRenders/BAPA/ColoredFragg-CBAPA_a-{attempt}_M-{M}_N-{N}_T-{Temp}.png")
#        render_scene(f"{data_directory}data/figures/aggRenders/BAPA/ColoredAsymFragg-{job_group}_a-{attempt}_N-{N}_T-{Temp}.png")
#        render_scene(f"{data_directory}data/figures/aggRenders/ColoredAgg-{job_group}_a-{attempt}_N-{N}.png")

    print('================================end================================')
    
    
    

#path = f'{data_directory}jobsCosine/lognorm_relax$a$/N_$n$/T_$t$/'
#attempt = 1
#N = 300
#T = 3
#path = path.replace('$a$',str(attempt)).replace('$n$',str(N)).replace('$t$',str(T))
paths = []
camera_dic = {}
color_dic = {}
energy_dic = {}


######################################################################################################################
#Lognorm and constant BPCA aggregates and camera dict
#paths.append(f'{data_directory}jobsCosine/lognorm_relax7/N_300/T_3/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax4/N_300/T_10/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax1/N_300/T_30/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax2/N_300/T_100/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax0/N_300/T_300/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax0/N_300/T_1000/')

#paths.append(f'{data_directory}jobsCosine/lognorm_relax6/N_100/T_3/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax1/N_100/T_10/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax7/N_100/T_30/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax0/N_100/T_100/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax5/N_100/T_300/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax3/N_100/T_1000/')

#paths.append(f'{data_directory}jobsCosine/lognorm_relax13/N_30/T_3/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax7/N_30/T_10/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax3/N_30/T_30/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax5/N_30/T_100/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax2/N_30/T_300/')
#paths.append(f'{data_directory}jobsCosine/lognorm_relax4/N_30/T_1000/')


#paths.append(f'{data_directory}jobsNovus/constrelax_1/N_300/T_3/')
#paths.append(f'{data_directory}jobsNovus/constrelax_9/N_300/T_10/')
#paths.append(f'{data_directory}jobsNovus/constrelax_17/N_300/T_30/')
#paths.append(f'{data_directory}jobsNovus/constrelax_5/N_300/T_100/')
#paths.append(f'{data_directory}jobsNovus/constrelax_4/N_300/T_300/')
#paths.append(f'{data_directory}jobsNovus/constrelax_1/N_300/T_1000/')

#paths.append(f'{data_directory}jobsNovus/constrelax_8/N_100/T_3/')
#paths.append(f'{data_directory}jobsNovus/constrelax_8/N_100/T_10/')
#paths.append(f'{data_directory}jobsNovus/constrelax_7/N_100/T_30/')
#paths.append(f'{data_directory}jobsNovus/constrelax_1/N_100/T_100/')
#paths.append(f'{data_directory}jobsNovus/constrelax_0/N_100/T_300/')
#paths.append(f'{data_directory}jobsNovus/constrelax_0/N_100/T_1000/')

#paths.append(f'{data_directory}jobsNovus/constrelax_18/N_30/T_3/')
#paths.append(f'{data_directory}jobsNovus/constrelax_18/N_30/T_10/')
#paths.append(f'{data_directory}jobsNovus/constrelax_5/N_30/T_30/')
#paths.append(f'{data_directory}jobsNovus/constrelax_3/N_30/T_100/')
#paths.append(f'{data_directory}jobsNovus/constrelax_7/N_30/T_300/')
#paths.append(f'{data_directory}jobsNovus/constrelax_2/N_30/T_1000/')



#camera_dic['/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_7/N_300/T_3/'] = (-1.6608,0.0000,0.6619) #0
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax4/N_300/T_10/'] = (1.2784,0.0000,-1.0415)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax1/N_300/T_30/'] = (-0.3552,0.0000,0.0127)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax2/N_300/T_100/'] = (2.0044,0.0000,1.5067)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax0/N_300/T_300/'] = (-1.5909,0.0000,-2.6332)
#camera_dic['/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_0/N_300/T_1000/'] = (-1.8004,0.0000,-0.2666)

#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax6/N_100/T_3/'] = (-1.6328,-3.1410,0.4030) #6
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax1/N_100/T_10/'] = (-1.3047,-0.0008,-0.0779)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax7/N_100/T_30/'] = (1.3831,3.1424,0.5578)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax0/N_100/T_100/'] = (0.9459,-0.2437,0.2140)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax5/N_100/T_300/'] = (-1.7864,-0.0008,1.3188)
#camera_dic['/mnt/be2a0173-321f-4b9d-b05a-addba547276f/SpaceLab_data/jobsCosine/lognorm_relax3/N_100/T_1000/'] = (-1.3117,-0.0008,-0.4339)

#camera_dic[f'{data_directory}jobsCosine/lognorm_relax13/N_30/T_3/'] = (1.6065,3.1424,1.2209) #6
#camera_dic[f'{data_directory}jobsCosine/lognorm_relax7/N_30/T_10/'] = (1.3971,3.1424,-0.8734)
#camera_dic[f'{data_directory}jobsCosine/lognorm_relax3/N_30/T_30/'] = (-1.6189,-0.0009,3.1478)
#camera_dic[f'{data_directory}jobsCosine/lognorm_relax5/N_30/T_100/'] = (-2.3170,-0.0007,-2.5066)
#camera_dic[f'{data_directory}jobsCosine/lognorm_relax2/N_30/T_300/'] = (1.2923,-0.0013,1.0202)
#camera_dic[f'{data_directory}jobsCosine/lognorm_relax4/N_30/T_1000/'] = (1.4878,-0.0013,-2.3940)


#camera_dic[f'{data_directory}jobsNovus/constrelax_1/N_300/T_3/'] = (1.7252,-0.0000,1.0878) #18
#camera_dic[f'{data_directory}jobsNovus/constrelax_9/N_300/T_10/'] = (0.6082,-0.0000,2.6517)
#camera_dic[f'{data_directory}jobsNovus/constrelax_17/N_300/T_30/'] = (-1.6747,-3.1416,-1.1811)
#camera_dic[f'{data_directory}jobsNovus/constrelax_5/N_300/T_100/'] = (-1.0045,-0.0000,0.7248)
#camera_dic[f'{data_directory}jobsNovus/constrelax_4/N_300/T_300/'] = (0.3708,-3.1416,-2.5425)
#camera_dic[f'{data_directory}jobsNovus/constrelax_1/N_300/T_1000/'] = (1.0131,-3.1416,0.5154)

#camera_dic[f'{data_directory}jobsNovus/constrelax_0/N_100/T_3/'] = (1.9556,-0.0001,0.1454) #24
#camera_dic[f'{data_directory}jobsNovus/constrelax_8/N_100/T_10/'] = (-1.2838,3.1417,-0.5387)
#camera_dic[f'{data_directory}jobsNovus/constrelax_7/N_100/T_30/'] = (1.3691,-0.0001,2.1561)
#camera_dic[f'{data_directory}jobsNovus/constrelax_1/N_100/T_100/'] = (1.2225,-0.0001,-2.1095)
#camera_dic[f'{data_directory}jobsNovus/constrelax_0/N_100/T_300/'] = (2.5001,-0.0001,1.3113)
#camera_dic[f'{data_directory}jobsNovus/constrelax_0/N_100/T_1000/'] = (-0.5784,0.1535,0.8012)

#camera_dic[f'{data_directory}jobsNovus/constrelax_18/N_30/T_3/'] = (-1.4024,3.1416,-0.0849) #30
#camera_dic[f'{data_directory}jobsNovus/constrelax_18/N_30/T_10/'] = (1.3552,-0.0001,1.9397)
#camera_dic[f'{data_directory}jobsNovus/constrelax_5/N_30/T_30/'] = (1.5925,-0.0001,-0.9506)
#camera_dic[f'{data_directory}jobsNovus/constrelax_3/N_30/T_100/'] =  (2.0184,-0.0000,1.5976)
#camera_dic[f'{data_directory}jobsNovus/constrelax_7/N_30/T_300/'] = (-0.9012,-0.0242,0.4328)
#camera_dic[f'{data_directory}jobsNovus/constrelax_2/N_30/T_1000/'] = (-0.3267,0.9415,-0.0829)
######################################################################################################################


######################################################################################################################
##Lognorm Asymmetric Fraggregates
#paths.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/AsymBAPA_3/N_300/T_1000/')
#camera_dic[paths[0]] = (1.5359, 0.0000, 2.4155)
##paths.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/AsymBAPA_4/N_300/T_1000/')
##camera_dic[paths[0]] = (-1.7904,-0.0041,0.9354)
#makeFrames(paths[-1],camera_dic[paths[-1]],render=True)
######################################################################################################################

######################################################################################################################
##BAPA for paper 2
dist = 100
emerald_green = '#008000'
energy = 5

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_1/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.3962665796279907,-3.0412348905883846e-07,-2.286388635635376]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_3/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.7904,-0.0041,0.9354]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_5/M_5/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.4668,-3.1416,0.0126]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_10/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.9051,1.4459e-7,0.9621]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_2/M_15/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-2.0029,8.489e-7,-2.0608]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_20/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.8074,6.7734e-7,1.7231]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_6/M_30/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.4862,9.9269e-7,1.2204]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_4/M_50/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.0604,1.2118e-6,-1.2789]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_60/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.8089,1.2434e-6,0.2919]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_75/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.1938034296035767,4.779301434609806e-07,2.3910927772521973]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_100/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-0.5282,-2.6997,-1.1036]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_150/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.6545690298080444,5.609205686596397e-07,-0.4433220624923706]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_300/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.3962665796279907,-3.0412348905883846e-07,-2.286388635635376]
# color_dic[path] = hex_to_rgb(emerald_green)
# energy_dic[path] = energy
#######################################################################################################################

######################################################################################################################
#CBAPA for paper 2
dist = 100
primary_blue = "#0000FF"
energy = 3

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_1/N_30/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist, 1.4451369047164917,-1.126496442793723e-07,-0.6527613401412964]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_3/N_90/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.6755133867263794,-1.4381688515641144e-06,-2.261948585510254]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_5/N_150/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.3962609767913818,-1.5987382084858837e-06,-2.5761075019836426]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_10/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.4800359010696411,-1.856592120930145e-06,-0.06981465965509415]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_15/N_450/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.3745945692062378,-0.7757140398025513,-0.4411630630493164]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_20/N_600/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-3.2578693208051845e-05,-1.4730552434921265,0.20947183668613434]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_1/M_30/N_900/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-0.0015141572803258896,-1.4800372123718262,1.6810728311538696]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_3/M_50/N_1500/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.3334301710128784,2.4294113245559856e-06,-3.019420862197876]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/CBAPA_3/M_60/N_1800/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.30551016330719,-1.9113298321826733e-07,-2.719231128692627]
# color_dic[path] = hex_to_rgb(primary_blue)
# energy_dic[path] = energy
######################################################################################################################

######################################################################################################################
#DBAPA for paper 2
dist = 100
electric_pink = "#FF1493"
energy = 3

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_20/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist, 1.3892877101898193,-4.4330184323371213e-07,2.349207878112793]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_40/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.4102312326431274,2.761186692623596e-07,-1.0157873630523682]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_60/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.4590998888015747,2.3489144496124936e-07,-1.3927781581878662]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_80/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.4660807847976685,5.77940966195456e-07,-1.49749755859375]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_100/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,0.607378363609314,4.1410379481021664e-07,0.1640557050704956]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_120/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-0.8726687431335449,3.6537102232614416e-07,2.5027971267700195]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_1/M_20/N_200/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.5847631692886353,8.123874550847177e-08,0.1919819861650467]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,2.4574272632598877,1.9517780458500056e-07,-0.4712454378604889]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_0/M_20/N_400/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.675519585609436,1.1326490323426697e-07,-0.205955371260643]
# color_dic[path] = hex_to_rgb(electric_pink)
# energy_dic[path] = energy

path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/DBAPA_9/M_20/N_600/T_1000/'
paths.append(path)
camera_dic[path] = [dist,0.18849888443946838,-1.3702323542474915e-07,1.4765418767929077]
color_dic[path] = hex_to_rgb(electric_pink)
energy_dic[path] = energy
######################################################################################################################

######################################################################################################################
#BAPAWELD for paper 2
dist = 100
yellowish_orange = "#FFA500"
energy = 2

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPAWELD_0/M_3/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,1.6755133867263794,-1.4381688515641144e-06,-2.261948585510254]
# color_dic[path] = hex_to_rgb(yellowish_orange)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPAWELD_3/M_15/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,0.5724671483039856,1.8504164245314314e-06,-1.0786194801330566]
# color_dic[path] = hex_to_rgb(yellowish_orange)
# energy_dic[path] = energy

# path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPAWELD_0/M_100/N_300/T_1000/'
# paths.append(path)
# camera_dic[path] = [dist,-1.528909683227539,2.9132195322745247e-06,1.5463529825210571]
# color_dic[path] = hex_to_rgb(yellowish_orange)
# energy_dic[path] = energy
######################################################################################################################

render = False
for key in camera_dic.keys():
    if len(camera_dic[key]) > 3:
        makeFrames(key,camera_dic[key][1:4],camera_dist=camera_dic[key][0],color=color_dic[key],sun_energy=energy_dic[key],render=render)
    elif len(camera_dic[key]) == 3:
        makeFrames(key,camera_dic[key],color=color_dic[key],sun_energy=energy_dic[key],render=render)

#    
    
    
#fragment figs
#makeFrames('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_20/N_300/T_1000/',camera_angle=(-0.2565,-0.9394,0.2276),render=True)
#makeFrames('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_1/M_30/N_300/T_1000/',camera_angle=(2.4953, -0.4932, -1.6288),render=True)
#makeFrames('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_50/N_300/T_1000/',camera_angle=(2.1048,-0.6182,0.3934),render=True)
#makeFrames('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_60/N_300/T_1000/',camera_angle=(2.3914,-0.5258,0.1591),render=True)
#makeFrames('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_100/N_300/T_1000/',camera_angle=(-0.5282,-2.6997,-1.1036),render=True)

#path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_20/N_300/T_1000/'
#makeFrames(path,render=True,camera_angle=(-1.9743,-0.5063,-0.2487))

#path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/SeqStickLognorm_20/N_300/'
#makeFrames(path,render=True,camera_angle=(-0.7077,-2.8008,0.8478))
#path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_7/N_300/T_3/'
#makeFrames(path,render=True,camera_angle=camera_dic[path])




#path = '/home/kolanzl/Desktop/Visualize/V2/'
#makeFrames(path,render=True,camera_angle=(2.1436,-0.4075,-1.9387))