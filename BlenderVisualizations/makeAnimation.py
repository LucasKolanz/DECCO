import bpy
import math
import os

def create_circle_animation(z_height=5, frames=100, total_time=10, radius=10, output_name='output'):
    fps = math.ceil((frames*1.0)/(total_time*1.0))
    print("fps: {}".format(fps))
    # Basic scene setup
    scene = bpy.context.scene
    scene.frame_end = frames
    scene.render.image_settings.file_format = 'FFMPEG'
    scene.render.fps = fps  # Set the frames per second here
    scene.render.ffmpeg.format = 'MPEG4'
    scene.render.ffmpeg.codec = 'H264'
    scene.render.ffmpeg.constant_rate_factor = 'LOW'#'MEDIUM'  # Quality setting

    # Ensure the output directory exists
    output_dir = bpy.path.abspath("//")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set the output file path for the animation
    scene.render.filepath = os.path.join(output_dir, f"{output_name}.mp4")

    # Check if there is a camera in the scene, if not, create one
    if not bpy.data.objects.get('Camera'):
        bpy.ops.object.camera_add(location=(0, 0, z_height))
    camera = bpy.data.objects['Camera']

    # Set the camera's display size
    camera.data.display_size = 0.5
    camera.data.clip_end = 10000


    # Create an empty object at the origin that the camera will track
    bpy.ops.object.empty_add(location=(0, 0, 0))
    origin_empty = bpy.data.objects['Empty']

    # Add a track to constraint to the camera so it always points at the empty object
    track_to_constraint = camera.constraints.new('TRACK_TO')
    track_to_constraint.target = origin_empty
    track_to_constraint.track_axis = 'TRACK_NEGATIVE_Z'
    track_to_constraint.up_axis = 'UP_Y'

    # Calculate the circular path for the camera
    for frame in range(frames):
        angle = 2*math.pi*frame/frames  # Current angle in radians
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)

        # Set camera position for current frame
        camera.location = (x, y, z_height)
        camera.keyframe_insert(data_path="location", frame=frame+1)

    # Render animation
    bpy.ops.render.render(animation=True)

    print(f"Animation created and saved as {output_name}.mp4")

# Call the function with desired parameters
create_circle_animation(z_height=500, frames=360, total_time=20, radius=1250, output_name='finalFigs/circleFrame11')
