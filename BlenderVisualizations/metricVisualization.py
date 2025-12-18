

import bpy
import sys
import os
import json
import tempfile
import numpy as np
from mathutils import Matrix



def addEnclosingEllipsoid(directory=""):
	# default to temp directory
	if directory == "":
		directory = tempfile.gettempdir()
	json_file = os.path.join(directory, "enclosing_ellipsoid.json")
	blender_ellipsoid_name = "EnclosingEllipsoid"

	# addEllipsoid(json_file,blender_ellipsoid_name,color=(0, 1, 0, 1))
	addEllipsoid(json_file,blender_ellipsoid_name,color=(255.0/256.0, 20.0/256.0, 147.0/256.0, 1))

def addEquivalentEllipsoid(directory=""):
	# default to temp directory
	if directory == "":
		directory = tempfile.gettempdir()
	json_file = os.path.join(directory, "equivalent_ellipsoid.json")
	blender_ellipsoid_name = "EquivalentEllipsoid"

	addEllipsoid(json_file,blender_ellipsoid_name,color=(0, 128.0/256.0 ,0,1))
	
def addEllipsoid(json_file, blender_ellipsoid_name="ellipsoid",color=(1,0,0,1)):
	import numpy as np
	import bpy
	from mathutils import Matrix

	print(f"Colro: {color}")

	scaleUp = 1e5  # same scaling used for spheres

	if not os.path.exists(json_file):
		print(f"[ERROR] file not found: {json_file}")
		return

	# --- Load ellipsoid data ---
	with open(json_file, "r") as f:
		data = json.load(f)

	center = np.array(data["center"], dtype=float)
	axes   = np.array(data["axes"],   dtype=float)  # (a,b,c)
	R      = np.array(data["R"],      dtype=float)  # (3x3)

	# Scale up for Blender
	center_bl = center * scaleUp
	axes_bl   = axes   * scaleUp

	print("[INFO] Loaded ellipsoid from JSON")
	print("       center (blend):", center_bl)
	print("       axes   (blend):", axes_bl)

	# --- Remove any existing object ---
	old = bpy.data.objects.get(blender_ellipsoid_name)
	if old:
		bpy.data.objects.remove(old, do_unlink=True)

	# --- Create a unit UV sphere at the origin ---
	bpy.ops.mesh.primitive_uv_sphere_add(
		radius=1.0,
		segments=64,
		ring_count=32,
		location=(0, 0, 0)
	)
	ellipsoid = bpy.context.active_object
	ellipsoid.name = blender_ellipsoid_name

	ellipsoid.scale = tuple(axes_bl)  # scale object data
	bpy.ops.object.transform_apply(scale=True)

	# Blender expects column-major rotation
	# --- Apply rotation and translation ---
	Rot = Matrix((
		(R[0,0], R[0,1], R[0,2]),
		(R[1,0], R[1,1], R[1,2]),
		(R[2,0], R[2,1], R[2,2])
	))

	# Convert to 4×4
	Rot4 = Rot.to_4x4()

	# Combine translation and rotation
	ellipsoid.matrix_world = Matrix.Translation(center_bl) @ Rot4


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
		bsdf.inputs["Base Color"].default_value = color

	# Assign material to mesh
	ellipsoid.data.materials.clear()
	ellipsoid.data.materials.append(mat)

	# Wireframe modifier settings to ensure material is used
	wire_mod = ellipsoid.modifiers.new(name="Wireframe", type='WIREFRAME')
	wire_mod.thickness = 0.1
	wire_mod.use_replace = True
	wire_mod.material_offset = 0

	# # --- Material ---
	# mat_name = "WireframeEllipsoidMaterial"
	# mat = bpy.data.materials.get(mat_name)
	# if mat is None:
	#     mat = bpy.data.materials.new(name=mat_name)

	# mat.use_nodes = True
	# bsdf = mat.node_tree.nodes.get("Principled BSDF")
	# if bsdf:
	#     bsdf.inputs["Base Color"].default_value = color

	# # CLEAR materials first
	# ellipsoid.data.materials.clear()

	# # Add a dummy base material (slot 0)
	# base_mat = bpy.data.materials.new(name="EllipsoidBasePlaceholder")
	# ellipsoid.data.materials.append(base_mat)

	# # Add your colored material (slot 1)
	# ellipsoid.data.materials.append(mat)

	# # Now wireframe modifier can reference slot 1
	# wf = ellipsoid.modifiers.new("Wireframe", type='WIREFRAME')
	# wf.thickness = 0.1
	# wf.use_replace = True
	# wf.material_offset = 1   # <-- use the colored one

	print("[INFO] Ellipsoid added successfully.")



def addEnclosingSphere(directory=""):
	# default to temp directory
	if directory == "":
		directory = tempfile.gettempdir()
	json_file = os.path.join(directory, "enclosing_sphere.json")
	blender_sphere_name = "EnclosingSphere"

	addSphere(json_file,blender_sphere_name,color=(191.0/256.0,0,191.0/256.0,1))

def addGyrationRadiusSphere(directory=""):
	# default to temp directory
	if directory == "":
		directory = tempfile.gettempdir()
	json_file = os.path.join(directory, "gyration_radius_sphere.json")
	blender_sphere_name = "GyrationRadiusSphere"

	addSphere(json_file,blender_sphere_name,color=(0,0,1,1))

def addSphere(json_file,blender_sphere_name="sphere",color=(1,0,0,1)):
	# Scale factor from simulation units to Blender units
	scaleUp = 1e5

	if not os.path.exists(json_file):
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
			bsdf.inputs["Base Color"].default_value = color

		# Assign material to mesh
		sphere.data.materials.clear()
		sphere.data.materials.append(mat)

		# Wireframe modifier settings to ensure material is used
		wire_mod = sphere.modifiers.new(name="Wireframe", type='WIREFRAME')
		wire_mod.thickness = 0.1
		wire_mod.use_replace = True
		wire_mod.material_offset = 0


def addConvexHull(directory="",blender_convexhull_name="convexHullMesh"):
	"""
	Loads the convex hull (written by write_convex_hull) and creates a
	Blender mesh object named 'ConvexHullMesh',or whatever is specified by blender_convexhull_name.

	The hull is visualized as a wireframe mesh.
	"""

	scaleUp = 1e5
	if directory == "":
		directory = tempfile.gettempdir()
	json_path = os.path.join(directory, "convex_hull.json")

	if not os.path.exists(json_path):
		print(f"[ERROR] Temp file not found: {json_path}")

	else:
		print(f"[INFO] Loading convex hull from: {json_path}")

		# --- Load JSON ---
		with open(json_path, "r") as f:
			data = json.load(f)

		# Option A: JSON contains "points"
		if "points" not in data:
			print("[ERROR] JSON does not contain 'points' field.")
			raise KeyError("Expected 'points' in convex_hull.json")

		points = np.array(data["points"], dtype=float)
		simplices = np.array(data["simplices"], dtype=int)

		# Scale points for Blender
		points_bl = points * scaleUp

		print(f"[INFO] Loaded {len(points_bl)} hull points")
		print(f"[INFO] Loaded {len(simplices)} hull faces")

		# --- Remove existing hull object if it exists ---
		HULL_NAME = blender_convexhull_name
		old = bpy.data.objects.get(HULL_NAME)
		if old is not None:
			bpy.data.objects.remove(old, do_unlink=True)

		# --- Create Blender mesh from points + triangles ---
		mesh = bpy.data.meshes.new(HULL_NAME)
		mesh.from_pydata(
			[tuple(p) for p in points_bl],
			[],
			[tuple(face) for face in simplices]
		)
		mesh.update()

		# Create Blender object
		obj = bpy.data.objects.new(HULL_NAME, mesh)
		bpy.context.collection.objects.link(obj)

		# --- Material ---
		mat_name = "WireframeHullMaterial"
		mat = bpy.data.materials.get(mat_name)
		if mat is None:
			mat = bpy.data.materials.new(name=mat_name)
			mat.use_nodes = True
			bsdf = mat.node_tree.nodes.get("Principled BSDF")
			bsdf.inputs["Base Color"].default_value = (191.0/256.0,191.0/256.0,0, 1.0)  # red

		obj.data.materials.clear()
		obj.data.materials.append(mat)

		# --- Wireframe modifier ---
		wire = obj.modifiers.new(name="Wireframe", type='WIREFRAME')
		wire.thickness = 0.1
		wire.use_replace = True
		wire.material_offset = 0  # <-- ensures the modifier uses our material



# ------------------------------------------------------------------------
# SUPER FAST shadow grid (single mesh object)
# ------------------------------------------------------------------------
def addShadowGrid(directory="", blender_shadow_grid_name="shadowGrid", grid_offset=10.0):
    """
    Super-fast version: builds ONE mesh with all *shadow* quads only.
    Cells where shadow[iy, ix] == False are simply not created.
    """
    scaleUp = 1e5
    if directory == "":
        directory = tempfile.gettempdir()
    json_path = os.path.join(directory, "shadow_grid.json")

    if not os.path.exists(json_path):
        print(f"[ERROR] File not found: {json_path}")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    xs     = np.asarray(data["xs"], float)
    ys     = np.asarray(data["ys"], float)
    shadow = np.asarray(data["shadow"], bool)
    delta  = float(data.get("delta", 0.01))
    k      = np.asarray(data.get("k", [0,0,1]), float)

    Nx = xs.shape[0]
    Ny = ys.shape[0]

    print(f"[INFO] Building fast shadow grid Nx={Nx}, Ny={Ny}")

    # Orthonormal basis
    u, v = make_uv_from_k(k)

    # Offset grid along +k
    offset_vec = k * (grid_offset)

    # Cell size in Blender units
    cell_size = delta * scaleUp
    h = cell_size / 2  # half size

    # ------------------------------------------------------------------
    # Build mesh data: ONLY for shadow == True
    # ------------------------------------------------------------------
    verts = []
    faces = []

    # plane corners in the (u, v) basis, centered at origin
    quad_local = [
        np.array([-h, -h]),
        np.array([ h, -h]),
        np.array([ h,  h]),
        np.array([-h,  h]),
    ]

    vid = 0
    for iy in range(Ny):
        y = ys[iy]
        for ix in range(Nx):
            # Skip non-shadow cells completely
            if not shadow[iy, ix]:
                continue

            x = xs[ix]

            # Cell center in plane coords
            P  = x*u + y*v
            P3 = P*scaleUp + offset_vec

            base = vid
            for q in quad_local:
                q3 = (P3
                      + q[0]*u
                      + q[1]*v)
                verts.append((q3[0], q3[1], q3[2]))
                vid += 1

            faces.append((base, base+1, base+2, base+3))

    print("[INFO] Vert count:", len(verts))
    print("[INFO] Face count:", len(faces))

    # ------------------------------------------------------------------
    # Create mesh object
    # ------------------------------------------------------------------
    mesh = bpy.data.meshes.new(blender_shadow_grid_name + "Mesh")
    obj  = bpy.data.objects.new(blender_shadow_grid_name, mesh)

    bpy.context.scene.collection.objects.link(obj)

    mesh.from_pydata(verts, [], faces)
    mesh.update()

    # ------------------------------------------------------------------
    # Single unlit material for all shadow quads
    # ------------------------------------------------------------------
    color = (140.0/256.0, 86.0/256.0, 75.0/256.0, 1.0)  # your brownish color
    grid_mat = make_simple_unlit_material("GridShadowUnlit", color)

    mesh.materials.clear()
    mesh.materials.append(grid_mat)

    # Make sure grid doesn't cast/receive shadows or affect lighting
    if hasattr(obj, "cycles_visibility"):
        vis = obj.cycles_visibility
        vis.shadow = False
        vis.diffuse = False
        vis.glossy = False
        vis.transmission = False
        vis.scatter = False

    if hasattr(obj, "hide_shadow"):
        obj.hide_shadow = True

    print("[INFO] Fast shadow grid created successfully!")


def make_simple_unlit_material(name, color=(1, 1, 1, 1)):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    # Emission-only shader (shadeless look)
    emission = nt.nodes.new("ShaderNodeEmission")
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    nt.links.new(emission.outputs["Emission"], out.inputs["Surface"])

    emission.inputs["Color"].default_value = color
    emission.inputs["Strength"].default_value = 1.0  # adjust if too bright/dim

    # Opaque, no alpha funny business
    mat.blend_method = "OPAQUE"
    mat.use_backface_culling = False

    return mat


# ------------------------------------------------------------------------
# Construct orthonormal basis (u, v) perpendicular to k
# ------------------------------------------------------------------------
def make_uv_from_k(k):
	k = np.asarray(k, float)
	k = k / np.linalg.norm(k)

	if abs(k[2]) < 0.9:
		a = np.array([0.0, 0.0, 1.0])
	else:
		a = np.array([0.0, 1.0, 0.0])

	u = np.cross(k, a)
	u /= np.linalg.norm(u)

	v = np.cross(k, u)
	return u, v
