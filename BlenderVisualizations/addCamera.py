import bpy
import mathutils

for camera in bpy.data.cameras:
    bpy.data.cameras.remove(camera)
        
bpy.ops.object.camera_add(location=(0, 0, 0))
    
# Set the newly created camera as the active camera
new_camera = bpy.context.object
bpy.context.scene.camera = new_camera

# Align the camera with the current view
bpy.ops.view3d.camera_to_view_selected()

# Set the camera's focal length to 1mm
new_camera.data.lens = 30.0


# Get the 3D view area
area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')

# Get the region data (this contains the view rotation)
region_data = next(space for space in area.spaces if space.type == 'VIEW_3D').region_3d

# Get the view direction
view_direction = region_data.view_rotation @ mathutils.Vector((0, 0, -1))

new_camera.rotation_euler = region_data.view_rotation.to_euler()

    # Move the camera backwards by 500 metersp
new_camera.location -= view_direction * 750
    
print("new_camera.rotation_euler: {}".format(new_camera.rotation_euler))
print("new_camera.location: {}".format(new_camera.location))