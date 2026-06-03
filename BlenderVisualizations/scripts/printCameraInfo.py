import bpy

#for i,camera in enumerate(bpy.data.cameras):
camera = bpy.data.objects['Camera']
i=0
#for i, camera in enumerate(bpy.data.objects['Camera']):
print("Camera {} info:".format(i))
print('angle   : ({},{},{})'.format(camera.rotation_euler.x,camera.rotation_euler.y,camera.rotation_euler.z))
print('location: {}'.format(camera.location))
print('focalLen: {}'.format(camera.data.lens))




