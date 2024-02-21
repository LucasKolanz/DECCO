import bpy

#for i,camera in enumerate(bpy.data.cameras):
camera = bpy.data.objects['Camera']
i=0
#for i, camera in enumerate(bpy.data.objects['Camera']):
print("Camera {} info:".format(i))
print('angle   : {}'.format(camera.rotation_euler))
print('location: {}'.format(camera.location))
print('focalLen: {}'.format(camera.data.lens))




