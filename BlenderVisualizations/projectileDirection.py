import bpy
import mathutils
import time
import fnmatch

scaleUp = 5e1

# Function to add an arrow in the direction of the given vector
def create_arrow(direction, name="ProjectileArrow"):
    # Create a new arrow mesh (cylinder and cone)
    bpy.ops.mesh.primitive_cylinder_add(radius=scaleUp*0.02, depth=scaleUp*0.5, location=(0, 0, 0))
    arrow_body = bpy.context.object
    arrow_body.name = name + "_body"

    # Move the body so the origin is at one end
    arrow_body.location = (0, 0, scaleUp*0.25)

    # Create cone for the arrow head
    bpy.ops.mesh.primitive_cone_add(radius1=scaleUp*0.05, depth=scaleUp*0.2, location=(0, 0, scaleUp*0.5))
    arrow_head = bpy.context.object
    arrow_head.name = name + "_head"

    # Combine arrow body and head into one object
    arrow_body.select_set(True)
    arrow_head.select_set(True)
    bpy.context.view_layer.objects.active = arrow_body
    bpy.ops.object.join()

    # Set the arrow's direction
    arrow_body.rotation_mode = 'QUATERNION'
    direction.normalize()
    arrow_body.rotation_quaternion = direction.to_track_quat('Z', 'Y')

    # Move the arrow to the origin
    arrow_body.location = (0, 0, 0)

# Function to extract projectile directions from a file
def read_projectile_directions(file_path):
    directions = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Projectile direction: "):
                # Extract the vector from the line (assuming it's formatted as "(x, y, z)")
                vector_str = line[len("Projectile direction: "):].strip()
                vector_str = vector_str.strip("()")
                vector_components = vector_str.split(",")
                direction = mathutils.Vector((float(vector_components[0]), 
                                              float(vector_components[1]), 
                                              float(vector_components[2])))
                directions.append(direction)
    return directions

# Main function to run in Blender
def visualize_projectile_directions(file_path):
    directions = read_projectile_directions(file_path)
    print(directions)
    for i, direction in enumerate(directions):
        time.sleep(1)
        create_arrow(direction, name=f"ProjectileArrow_{i}")
        
def visualize_projectile_directions_with_delay(file_path, delay=0.5):
    directions = read_projectile_directions(file_path)
    
    def show_arrow(index=0):
        if index < len(directions):
            create_arrow(directions[index], name=f"ProjectileArrow_{index}")
            # Schedule the next arrow after the delay, incrementing the index
            bpy.app.timers.register(lambda: show_arrow(index + 1), first_interval=delay)
        return None  # Stop the timer when done
    
    # Start by showing the first arrow
    show_arrow(0)
    
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
    
# Specify the file path
#file_path = "/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA0/N_300/T_1000/sim_errors.txt"
#file_path = "/media/kolanzl/easystore/SpaceLab_data/jobs/TESTBAPA1/N_300/T_1000/sim_errors.txt"
#file_path = '/media/kolanzl/easystore/SpaceLab_data/jobs/TESTBAPA2/N_3000/T_1000/sim_errors.txt'
file_path = '/media/kolanzl/easystore/SpaceLab_data/jobs/TESTBAPA3/N_3000/T_1000/sim_errors.txt'
file_path = '/media/kolanzl/easystore/SpaceLab_data/jobs/TESTBPCA0/N_30/T_3/sim_errors.txt'
file_path = '/media/kolanzl/easystore/SpaceLab_data/jobs/BAPA8/N_300/T_1000/sim_errors.txt'


# Visualize the projectile directions
visualize_projectile_directions_with_delay(file_path)
