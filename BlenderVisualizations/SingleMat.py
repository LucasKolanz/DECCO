from __future__ import division
import bpy
import numpy as np
#from mathutils import *
from math import *
import fnmatch
import os
import linecache
import mathutils

#def create_animation(frames=100, total_time=10, output_name='output'):
def create_animation(output_name='output',frames=10):
#    fps = math.ceil((frames*1.0)/(total_time*1.0))
#    print("fps: {}".format(fps))
    # Basic scene setup
    scene = bpy.context.scene
    scene.frame_end = frames
    scene.render.image_settings.file_format = 'FFMPEG'
#    scene.render.fps = fps  # Set the frames per second here
    scene.render.ffmpeg.format = 'MPEG4'
    scene.render.ffmpeg.codec = 'H264'
    scene.render.ffmpeg.constant_rate_factor = 'LOW'#'MEDIUM'  # Quality setting

    # Ensure the output directory exists
    output_dir = bpy.path.abspath("//")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set the output file path for the animation
    scene.render.filepath = os.path.join(output_dir, f"{output_name}.mp4")

    # Render animation
    bpy.ops.render.render(animation=True)

    print(f"Animation created and saved as {output_name}.mp4")


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
    render = bpy.context.scene.render
    render.image_settings.file_format = 'PNG'  # or 'JPEG', 'BMP', etc.
    render.filepath = picName  # Change this to your desired path
    
    # Render the image
    bpy.ops.render.render(write_still=True)

def add_camera(camera_angle='default',camera_position='default',shift=(0,0),focal_length=39.0):
#    delete_cameras()

    print("Camera info in add_camera: {}, {}".format(camera_angle,camera_position))
    
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
    
    
    if camera_angle == 'default' and camera_position == 'default':
        # Move the camera back along the view vector
        #        mat = new_camera.matrix_world

        # Get the 3D view area
        area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')

        # Get the region data (this contains the view rotation)
        region_data = next(space for space in area.spaces if space.type == 'VIEW_3D').region_3d

        # Get the view direction
        view_direction = region_data.view_rotation @ mathutils.Vector((0, 0, -1))

        new_camera.rotation_euler = region_data.view_rotation.to_euler()

            # Move the camera backwards by 500 metersp
        new_camera.location -= view_direction * 750
    else:
        new_camera.rotation_euler = camera_angle
            # Move the camera backwards by 500 metersp
        new_camera.location = camera_position

    
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
    principled_bsdf.inputs['Base Color'].default_value = (0.01, 0.01,0.01, opacity)  # RGBA values, all set to 0 for black

    # Add Material Output node
    material_output = shader_tree.nodes.new(type='ShaderNodeOutputMaterial')
    material_output.location = (300, 0)
    shader_tree.links.new(principled_bsdf.outputs['BSDF'], material_output.inputs['Surface'])
    
    material.blend_method = 'HASHED'
    material.shadow_method = 'HASHED'
    material.use_screen_refraction = True

    return material

def NthOpacity(N,opacity):
    index = np.argpartition(opacity,-N)[-N]
    return opacity[index]


def makeFrames(frames,camera_angle = "default",camera_position = "default",camera_name = "default",shift=(0,0),focal_length=39.0):
    if type(frames) != list:
        frames = [frames]
        
    frameNum = 0
    
    
        
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

    set_background()
    add_sun()
    if camera_angle != "default" and camera_position != "default":
#        camera = ''
        print('deleting cameras')
        delete_cameras()
#        print("Camera info from doEverything: {}, {}, {}".format(camera_angle,camera_position,camera_name))
        camera=add_camera(camera_angle,camera_position,shift,focal_length)


             
    path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/JobPaperData/'
    simFile = path + "1250_5000-n6.0-IP22-dt0.04_simData.csv"
    constFile = path + "1250_5000-n9.0-IP22-dt0.04_constants.csv"
#    opacityFile = path + "OpacityValuesT.txt"
    opacityFile = path + "allStrainT.txt"

    properties = 11
    scaleUp = 1e-5
    
    lowOpacity = 0.03
    highOpacity = 0.55
    sphereSet = []

    
    for frame in frames:
        print('===============================start frame{}==============================='.format(frame))



        
        #input()
    #    print("HEREREREREREREREEERER: {}".format(camera_angle != "default" and camera_position != "default"))
       
        
        
        

       

        print(frame)


#        print(linecache.getline(simFile,frame).strip(' \n').split(','))
#        print(linecache.getline(opacityFile,frame+1).strip(' \n').split(' '))
        simData = np.array(linecache.getline(simFile,frame+2).strip(' \n').split(','),dtype=np.float64)
#        opacityData1 = np.array(linecache.getline(opacityFile,frame+1).strip(' \n').split(' '),dtype=np.float64)
        opacityData = np.loadtxt(opacityFile,delimiter = ' ')[frame,:]#np.array(linecache.getline(path+'allStrainT.txt',1).strip(' \n').split(' '),dtype=np.float64)
#        print(opacityData1[0:5])
#        print(opacityData[0:5])
#        print("///////////////////////////////////////////////////////////////")
#        print(np.array(linecache.getline(opacityFile,frame+1).strip('\n').split(' ')).shape)
#        print(opacityData.shape)
        opacityMax = opacityData[opacityData!=0].max()
        opacityStd = opacityData[opacityData!=0].std()
        opacityCutoff = NthOpacity(75,opacityData)
        
        print("HERERE")
        print(opacityCutoff)

        opacityData = np.where(opacityData < opacityCutoff,0.0,1.0)
        print(f"There are {np.sum(opacityData)} high opacity values in frame {frame}")
        #opacityData = np.where(opacityData > opacityMean + opacityStd,0,opacityData)
        #opacityData = np.where(opacityData > 0,1,0)

        constData = np.loadtxt(constFile,dtype=float,delimiter=',',skiprows = 0)

        #####################################################
        #constData = constData[opacityData>0,:]
        #temp_simData = np.zeros(constData.shape[0]*11)
        #print(temp_simData.shape)

        #i = 0
        #for ind,data in enumerate(opacityData):
        #    if data != 0:
        #        temp_simData[i*11:(i+1)*11] = simData[ind*11:(ind+1)*11]
        #        i += 1
        #    
        #simData = temp_simData
        #opacityData = opacityData[opacityData != 0]
        #####################################################

        print("simdata shape: " + str(simData.shape))
        print("opacityData shape: " + str(opacityData.shape))
        print("constData shape: " + str(constData.shape))

        numSpheres = constData.shape[0]

        num=numSpheres
        numSpheres = num
        simData = simData[:num*11]
        opacityData = opacityData[:num]
        constData = constData[:num,:]


        #print("simdata shape: " + str(simData.shape))
        
        ## Initial sphere mesh to be instanced:
        print('create mesh')
        bpy.ops.mesh.primitive_uv_sphere_add(location = (0,0,0), radius = 1)
        obj = bpy.context.object # the currently selected object
        sphereMesh = obj.data # retrieve the mesh
        bpy.data.objects.remove(obj) # remove the object


        print("num spheres: ",numSpheres)
        


        #make material
    #    materials = []
        transparentMaterial = create_black_material(0)
        transparentMaterial.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = lowOpacity
        opaqueMaterial = create_black_material(1)
        opaqueMaterial.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = highOpacity


        print('Instanciate spheres')
        # Instanciate spheres:
#        bpy.context.scene.frame_set(frameNum)
        if len(sphereSet) == 0:
            for sphere_num in range(numSpheres):
                sphere = bpy.data.objects.new("Mball." + str(sphere_num),sphereMesh)
                if not sphere.material_slots:
            #        with bpy.context.temp_override({'object': sphere}):
                    bpy.ops.object.material_slot_add({'object': sphere})
            #    print("opacity: "+str(opacityData[sphere_num]))
        #        materials.append(create_black_material(sphere_num))
                sphere.material_slots[0].link = 'OBJECT'
                bpy.context.scene.collection.objects.link(sphere) # link the object to the scene collection
                sphereSet.append(sphere)

        for sphere_num in range(numSpheres):
#            sphere = sphereSet[sphere_num]

            
            if opacityData[sphere_num] == 1:
                mat = opaqueMaterial
            else:
                mat = transparentMaterial
            sphereSet[sphere_num].material_slots[0].material = mat
    #        materials[sphere_num].node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = opacityData[sphere_num]
            sphereSet[sphere_num].scale = (scaleUp*constData[sphere_num,0],scaleUp*constData[sphere_num,0],scaleUp*constData[sphere_num,0])
            sphereSet[sphere_num].location = (scaleUp*simData[properties*sphere_num+0],scaleUp*simData[properties*sphere_num+1],scaleUp*simData[properties*sphere_num+2])
            sphereSet[sphere_num].rotation_mode = "XYZ"
            
        if len(frames) > 1:
            camera.keyframe_insert(data_path="location", frame=-1)
                
            frameNum += 1                
            
#            
#    if len(frames) == 1:
        render_scene("finalFigs/SingleShots/frame_{:04d}-view_{}.png".format(frames[0],camera_name))
#    elif len(frames) > 1:
#        create_animation(output_name="finalFigs/animation_{}-view_{}.png".format(frames,camera_name),frames=frameNum)
    print('================================end================================')
    
    
    
camera_angles = [(-3.0103, -0.0013, 0.7735),(1.2357, 0.0000, 1.4591)]
camera_locations = [(-90.5415, 94.5668, -991.3929),(938.4911, -105.2690, 328.8664)]
camera_names = ['stringy','side']

stillFrame25StringyCam = [(-3.0103, -0.0013, 0.7735),(-90.5415, 94.5668, -991.3929),'stringy',(0.04,0.0),30]#[angle,location,name,shift,focal_length]
    
stillFrame0sideCam = [(1.2357, 0.0000, 1.4591),(938.4911, -105.2690, 328.8664),'side',(0.04,-0.04),36]

li = [0,2,10,25]
li = list(range(0,25))
#li = [2]

#makeFrames(25,stillFrame25StringyCam[0],stillFrame25StringyCam[1],stillFrame25StringyCam[2],stillFrame25StringyCam[3],stillFrame25StringyCam[4])
#makeFrames(25,camera_angles[0],camera_locations[0],camera_names[0])
makeFrames(1,stillFrame0sideCam[0],stillFrame0sideCam[1],'camtest_side',stillFrame0sideCam[3],stillFrame0sideCam[4])

#for c in range(len(camera_angles)):
#    for frame in li:
#        makeFrames(frame,camera_angles[c],camera_locations[c],camera_names[c],(0.04,0))


#makeFrames(list(range(25)),camera_angles[c],camera_locations[c],camera_names[c])
#makeFrames([0,1],camera_angles[c],camera_locations[c],camera_names[c])