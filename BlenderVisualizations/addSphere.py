


import bpy
import json
import os
import tempfile

def addEnclosingSphere(directory=""):
    # default to temp directory
    if directory == "":
        directory = tempfile.gettempdir()
    json_file = os.path.join(directory, "enclosing_sphere.json")
    blender_sphere_name = "EnclosingSphere"

    addSphere(json_file,blender_sphere_name)

def addGyrationRadiusSphere(directory=""):
    # default to temp directory
    if directory == "":
        directory = tempfile.gettempdir()
    json_file = os.path.join(directory, "gyration_radius_sphere.json")
    blender_sphere_name = "GyrationRadiusSphere"

    addSphere(json_file,blender_sphere_name)

def addSphere(json_file,blender_sphere_name="sphere"):
    # Scale factor from simulation units to Blender units
    scaleUp = 1e5

    if not os.path.exists(json_path):
        print(f"[ERROR] File not found: {json_file}")
    else:
        # --- Read center & radius from JSON ---
        with open(json_file, "r") as f:
            data = json.load(f)

        center = data["center"]   # [x, y, z] in simulation units
        r      = data["radius"]   # scalar in simulation units

        # Apply the same scaleUp you’ve been using
        center = [c * scaleUp for c in center]
        r = r * scaleUp

        print(f"[INFO] Loaded sphere from {json_file}")
        print(f"       center = {center}")
        print(f"       radius = {r}")

        # Optional: clean up existing object named "BoundingSphere"
        # so you don’t stack multiple ones.
        # old_obj = bpy.data.objects.get("Sphere")
        # if old_obj is not None:
        #     bpy.data.objects.remove(old_obj, do_unlink=True)

        # --- Add a UV sphere ---
        bpy.ops.mesh.primitive_uv_sphere_add(
            radius=r,
            segments=64,
            ring_count=32,
            location=center
        )

        sphere = bpy.context.active_object
        sphere.name = "BoundingSphere"

        # Smooth shading (optional)
        bpy.ops.object.shade_smooth()

        # --- Create and assign a material (optional) ---
        mat_name = "WireframeMaterial"
        mat = bpy.data.materials.get(mat_name)
        if mat is None:
            mat = bpy.data.materials.new(name=mat_name)
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf:
            bsdf.inputs["Base Color"].default_value = (1.0, 0.0, 0.0, 1.0)

        # Assign material to mesh
        ellipsoid.data.materials.clear()
        ellipsoid.data.materials.append(mat)

        # Wireframe modifier settings to ensure material is used
        wire_mod = ellipsoid.modifiers.new(name="Wireframe", type='WIREFRAME')
        wire_mod.thickness = 0.01
        wire_mod.use_replace = True
        wire_mod.material_offset = 0
