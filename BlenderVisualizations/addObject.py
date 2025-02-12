import bpy
import sys
project_path = '/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab/'
sys.path.append(project_path+"utilities/")
sys.path.append('/home/kolanzl/.local/lib/python3.11/site-packages')
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u
import numpy as np

def calc_porosity_KBM(data_folder,size):
    data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=False)
    if data is None:
        return np.nan
    # num_balls = data.shape[0]

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3)  
    # effective_radius = radius*np.power(num_balls,1/3) 
        
    principal_moi = u.get_principal_moi(np.mean(mass),data)
    # principal_moi = get_principal_moi(mass,data)

    alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)
    # alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)

    RKBM = np.sqrt(np.sum(alphai)/3) * effective_radius
    return RKBM

def calc_porosity_abc(data_folder,size):
    data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=False)
    if data is None:
        return np.nan
    # num_balls = data.shape[0]

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 


    # effective_radius = radius*np.power(num_balls,1/3) 
        
    principal_moi = u.get_principal_moi(np.mean(mass),data)
    # principal_moi = get_principal_moi(mass,data)
    
    
    alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)
    # alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)
    
    a = effective_radius * np.sqrt(alphai[1] + alphai[2] - alphai[0])
    b = effective_radius * np.sqrt(alphai[2] + alphai[0] - alphai[1])
    c = effective_radius * np.sqrt(alphai[0] + alphai[1] - alphai[2])
    
    # Rabc = np.power(a*b*c,1/3)
#    porosity = 1-(effective_radius**3/(a*b*c))
    return a,b,c

def addCircle(path,num_particles):
    scaleUp = 1e5
    r = calc_porosity_KBM(path,num_particles)*scaleUp
    ## Clean up existing objects (optional)
    #bpy.ops.object.select_all(action='SELECT')
    #bpy.ops.object.delete()

    # Add a UV sphere
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=1.0,
        segments=32,
        ring_count=16,
        location=(0,0,0)
    )

    sphere = bpy.context.active_object
    sphere.scale = (r,r,r)

    # Smooth shading (optional)
    bpy.ops.object.shade_smooth()

    # Add a wireframe modifier
    bpy.ops.object.modifier_add(type='WIREFRAME')
    wire_mod = sphere.modifiers["Wireframe"]

    # Adjust the thickness of the wire
    wire_mod.thickness = 0.01

    # (Optional) Apply the modifier if you want the mesh to be permanently wireframed
    # bpy.ops.object.modifier_apply(modifier="Wireframe")

    # Create and assign a material (optional)
    mat = bpy.data.materials.new(name="circleWireframeMaterial")
    mat.diffuse_color = (0.0, 0.0, 1.0, 1)  #  color
    sphere.data.materials.append(mat)


def addEllipse(path,num_particles):
    # Define your semi-major axes
    # For example, if your ellipsoid has semi-major axes a, b, and c along x, y, z respectively:
    scaleUp = 1e5
    a,b,c=calc_porosity_abc(path,num_particles)
    a=a*scaleUp
    b=b*scaleUp
    c=c*scaleUp

    ## Clean up existing objects (optional)
    #bpy.ops.object.select_all(action='SELECT')
    #bpy.ops.object.delete()

    # Add a UV sphere
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=1.0,
        segments=32,
        ring_count=16,
        location=(0,0,0)
    )

    ellipsoid = bpy.context.active_object
    ellipsoid.scale = (a, b, c)

    # Smooth shading (optional)
    bpy.ops.object.shade_smooth()

    # Add a wireframe modifier
    bpy.ops.object.modifier_add(type='WIREFRAME')
    wire_mod = ellipsoid.modifiers["Wireframe"]

    # Adjust the thickness of the wire
    wire_mod.thickness = 0.01

    # (Optional) Apply the modifier if you want the mesh to be permanently wireframed
    # bpy.ops.object.modifier_apply(modifier="Wireframe")

    # Create and assign a material (optional)
    mat = bpy.data.materials.new(name="WireframeMaterial")
    mat.diffuse_color = (1.0, 0.0, 0.0, 1)  #  color
    ellipsoid.data.materials.append(mat)
