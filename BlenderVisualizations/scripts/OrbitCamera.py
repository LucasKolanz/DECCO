#Run this after loading in a simulation and you can generate
#an image sequence of the simulation to stitch together with
#FFMPEG



import bpy
import math
from mathutils import Vector
import os
import random

# def sphere_project_uv_even_if_hidden(obj):
def sphere_project_uv(obj):
    """
    Create spherical UV projection for one mesh object.
    Temporarily unhides the object if needed, then restores visibility.
    """

    if obj.type != "MESH":
        return False

    if obj.name not in bpy.context.view_layer.objects:
        print(f"Skipping {obj.name}: not in current view layer")
        return False

    old_hide_viewport = obj.hide_viewport
    old_hide_render = obj.hide_render
    old_hide_get = obj.hide_get()

    try:
        # Temporarily make object editable
        obj.hide_viewport = False
        obj.hide_render = False
        obj.hide_set(False)

        bpy.ops.object.select_all(action='DESELECT')
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.sphere_project()
        bpy.ops.object.mode_set(mode='OBJECT')

        return True

    except RuntimeError as err:
        print(f"UV projection failed for {obj.name}: {err}")
        try:
            bpy.ops.object.mode_set(mode='OBJECT')
        except RuntimeError:
            pass
        return False

    finally:
        # Restore animation/render visibility state
        obj.hide_viewport = old_hide_viewport
        obj.hide_render = old_hide_render
        obj.hide_set(old_hide_get)
        obj.select_set(False)

def set_object_material(obj, material):
    """
    Assign a material to a Blender object.
    """
    if obj.data.materials:
        obj.data.materials[0] = material
    else:
        obj.data.materials.append(material)

import bpy
import os
import math

def make_shared_random_oriented_equirectangular_material(
    image_path,
    name="SharedRandomEquirectMaterial",
    roughness=0.8,
    metallic=0.0,
    specular=0.5,
    max_rotation_degrees=360.0,
):
    """
    Create ONE shared material that uses Object Info -> Random
    to rotate the equirectangular texture differently on each object.

    This assumes the spheres have UVs.
    """

    if not os.path.isfile(image_path):
        raise FileNotFoundError(f"Image not found: {image_path}")

    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)

    mat.use_nodes = True

    nt = mat.node_tree
    nodes = nt.nodes
    links = nt.links
    nodes.clear()

    # Nodes
    texcoord = nodes.new(type="ShaderNodeTexCoord")
    texcoord.location = (-1200, 0)

    obj_info = nodes.new(type="ShaderNodeObjectInfo")
    obj_info.location = (-1200, -250)

    # Move UV center to origin so rotation happens around center of image
    vec_sub = nodes.new(type="ShaderNodeVectorMath")
    vec_sub.operation = 'SUBTRACT'
    vec_sub.location = (-950, 0)
    vec_sub.inputs[1].default_value = (0.5, 0.5, 0.0)

    # Convert random number in [0,1] to angle in radians
    angle_mult = nodes.new(type="ShaderNodeMath")
    angle_mult.operation = 'MULTIPLY'
    angle_mult.location = (-950, -250)
    angle_mult.inputs[1].default_value = math.radians(max_rotation_degrees)

    # Rotate UV vector around Z
    vec_rotate = nodes.new(type="ShaderNodeVectorRotate")
    vec_rotate.location = (-700, 0)
    vec_rotate.rotation_type = 'AXIS_ANGLE'
    vec_rotate.inputs["Axis"].default_value = (0.0, 0.0, 1.0)
    vec_rotate.inputs["Center"].default_value = (0.0, 0.0, 0.0)

    # Move UV center back
    vec_add = nodes.new(type="ShaderNodeVectorMath")
    vec_add.operation = 'ADD'
    vec_add.location = (-450, 0)
    vec_add.inputs[1].default_value = (0.5, 0.5, 0.0)

    image_tex = nodes.new(type="ShaderNodeTexImage")
    image_tex.location = (-200, 0)
    image_tex.image = bpy.data.images.load(image_path, check_existing=True)
    image_tex.image.colorspace_settings.name = 'sRGB'

    bsdf = nodes.new(type="ShaderNodeBsdfPrincipled")
    bsdf.location = (50, 0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic

    if "Specular IOR Level" in bsdf.inputs:
        bsdf.inputs["Specular IOR Level"].default_value = specular
    elif "Specular" in bsdf.inputs:
        bsdf.inputs["Specular"].default_value = specular

    out = nodes.new(type="ShaderNodeOutputMaterial")
    out.location = (300, 0)

    # Links
    links.new(texcoord.outputs["UV"], vec_sub.inputs[0])
    links.new(obj_info.outputs["Random"], angle_mult.inputs[0])

    links.new(vec_sub.outputs["Vector"], vec_rotate.inputs["Vector"])
    links.new(angle_mult.outputs["Value"], vec_rotate.inputs["Angle"])

    links.new(vec_rotate.outputs["Vector"], vec_add.inputs[0])
    links.new(vec_add.outputs["Vector"], image_tex.inputs["Vector"])

    links.new(image_tex.outputs["Color"], bsdf.inputs["Base Color"])
    links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    return mat

def object_has_uvs(obj):
    return obj.type == "MESH" and len(obj.data.uv_layers) > 0

def apply_equirectangular_texture_to_spheres(
    image_path,
    material_name="DustEquirectangularMaterial",
    sphere_name_patterns=("Sphere", "sphere", "Ball", "ball"),
    roughness=0.8,
    metallic=0.0,
    do_uv_project=True,
):
    """
    Apply an equirectangular image texture to all sphere-like mesh objects.
    """

    mat = make_equirectangular_sphere_material(
        image_path=image_path,
        name=material_name,
        roughness=roughness,
        metallic=metallic,
    )

    count = 0

    for obj in bpy.context.scene.objects:
        if obj.type != "MESH":
            continue

        if any(pattern in obj.name for pattern in sphere_name_patterns):
            if do_uv_project and not object_has_uvs(obj):
                sphere_project_uv(obj)

            set_object_material(obj, mat)
            count += 1

    print(f"Applied equirectangular texture to {count} sphere-like objects.")

    return mat

def clear_world_background(color=(0.0, 0.0, 0.0, 1.0), strength=1.0):
    """
    Reset the world background to a plain color.
    This clears any previous environment texture setup.
    """

    scene = bpy.context.scene

    if scene.world is None:
        scene.world = bpy.data.worlds.new("World")

    world = scene.world
    world.use_nodes = True

    nt = world.node_tree
    nodes = nt.nodes
    links = nt.links

    nodes.clear()

    out = nodes.new(type="ShaderNodeOutputWorld")
    out.location = (300, 0)

    bg = nodes.new(type="ShaderNodeBackground")
    bg.location = (0, 0)
    bg.inputs["Color"].default_value = color
    bg.inputs["Strength"].default_value = strength

    links.new(bg.outputs["Background"], out.inputs["Surface"])

    scene.render.film_transparent = False


def set_equirectangular_background(
    image_path,
    strength=1.0,
    rotation_x_deg=0.0,
    rotation_y_deg=0.0,
    rotation_z_deg=0.0,
):
    """
    Use an equirectangular image as the world background.

    image_path:
        Absolute path to an equirectangular image file.

    strength:
        Background/light intensity.

    rotation_*_deg:
        Rotate the environment map in degrees.
        Usually rotation_z_deg is the main one you want to change.
    """

    if not os.path.isfile(image_path):
        raise FileNotFoundError(f"Background image not found: {image_path}")

    scene = bpy.context.scene

    if scene.world is None:
        scene.world = bpy.data.worlds.new("World")

    world = scene.world
    world.use_nodes = True

    nt = world.node_tree
    nodes = nt.nodes
    links = nt.links

    nodes.clear()

    texcoord = nodes.new(type="ShaderNodeTexCoord")
    texcoord.name = "BG_TexCoord"
    texcoord.location = (-900, 0)

    mapping = nodes.new(type="ShaderNodeMapping")
    mapping.name = "BG_Mapping"
    mapping.location = (-700, 0)
    mapping.inputs["Rotation"].default_value[0] = math.radians(rotation_x_deg)
    mapping.inputs["Rotation"].default_value[1] = math.radians(rotation_y_deg)
    mapping.inputs["Rotation"].default_value[2] = math.radians(rotation_z_deg)

    env = nodes.new(type="ShaderNodeTexEnvironment")
    env.name = "BG_Environment"
    env.location = (-450, 0)
    env.image = bpy.data.images.load(image_path, check_existing=True)
    env.image.colorspace_settings.name = 'sRGB'

    bg = nodes.new(type="ShaderNodeBackground")
    bg.name = "BG_Background"
    bg.location = (-150, 0)
    bg.inputs["Strength"].default_value = strength

    out = nodes.new(type="ShaderNodeOutputWorld")
    out.name = "BG_Output"
    out.location = (150, 0)

    links.new(texcoord.outputs["Generated"], mapping.inputs["Vector"])
    links.new(mapping.outputs["Vector"], env.inputs["Vector"])
    links.new(env.outputs["Color"], bg.inputs["Color"])
    links.new(bg.outputs["Background"], out.inputs["Surface"])

    scene.render.film_transparent = False


def set_equirectangular_background_visible_only(
    image_path,
    camera_strength=1.0,
    lighting_strength=0.0,
    rotation_x_deg=0.0,
    rotation_y_deg=0.0,
    rotation_z_deg=0.0,
):
    """
    Use an equirectangular image as the visible background, while controlling
    how much it lights the scene.

    camera_strength:
        How bright the image looks in the background.

    lighting_strength:
        How much light it contributes to the scene.
        Use 0.0 for 'background visible but does not light the scene'.
    """

    if not os.path.isfile(image_path):
        raise FileNotFoundError(f"Background image not found: {image_path}")

    scene = bpy.context.scene

    if scene.world is None:
        scene.world = bpy.data.worlds.new("World")

    world = scene.world
    world.use_nodes = True

    nt = world.node_tree
    nodes = nt.nodes
    links = nt.links

    nodes.clear()

    texcoord = nodes.new(type="ShaderNodeTexCoord")
    texcoord.location = (-1100, 0)

    mapping = nodes.new(type="ShaderNodeMapping")
    mapping.location = (-900, 0)
    mapping.inputs["Rotation"].default_value[0] = math.radians(rotation_x_deg)
    mapping.inputs["Rotation"].default_value[1] = math.radians(rotation_y_deg)
    mapping.inputs["Rotation"].default_value[2] = math.radians(rotation_z_deg)

    env = nodes.new(type="ShaderNodeTexEnvironment")
    env.location = (-650, 0)
    env.image = bpy.data.images.load(image_path, check_existing=True)
    env.image.colorspace_settings.name = 'sRGB'

    bg_camera = nodes.new(type="ShaderNodeBackground")
    bg_camera.location = (-350, 120)
    bg_camera.inputs["Strength"].default_value = camera_strength

    bg_lighting = nodes.new(type="ShaderNodeBackground")
    bg_lighting.location = (-350, -120)
    bg_lighting.inputs["Strength"].default_value = lighting_strength

    light_path = nodes.new(type="ShaderNodeLightPath")
    light_path.location = (-350, -320)

    mix = nodes.new(type="ShaderNodeMixShader")
    mix.location = (-80, 0)

    out = nodes.new(type="ShaderNodeOutputWorld")
    out.location = (150, 0)

    links.new(texcoord.outputs["Generated"], mapping.inputs["Vector"])
    links.new(mapping.outputs["Vector"], env.inputs["Vector"])

    links.new(env.outputs["Color"], bg_camera.inputs["Color"])
    links.new(env.outputs["Color"], bg_lighting.inputs["Color"])

    # 0 -> lighting background, 1 -> camera background
    links.new(light_path.outputs["Is Camera Ray"], mix.inputs["Fac"])
    links.new(bg_lighting.outputs["Background"], mix.inputs[1])
    links.new(bg_camera.outputs["Background"], mix.inputs[2])

    links.new(mix.outputs["Shader"], out.inputs["Surface"])

    scene.render.film_transparent = False


def rotate_equirectangular_background(rotation_z_deg):
    """
    Rotate an already-created world environment around Z.
    Useful for trying different sky orientations quickly.
    """

    world = bpy.context.scene.world
    if world is None or not world.use_nodes:
        raise RuntimeError("Scene world is not using nodes.")

    mapping = world.node_tree.nodes.get("BG_Mapping")
    if mapping is None:
        raise RuntimeError("Could not find BG_Mapping node.")

    mapping.inputs["Rotation"].default_value[2] = math.radians(rotation_z_deg)


def set_world_background_strength(strength):
    """
    Change the strength of a simple equirectangular world background.
    """

    world = bpy.context.scene.world
    if world is None or not world.use_nodes:
        raise RuntimeError("Scene world is not using nodes.")

    bg = world.node_tree.nodes.get("BG_Background")
    if bg is None:
        raise RuntimeError("Could not find BG_Background node.")

    bg.inputs["Strength"].default_value = strength


def clear_equirectangular_background():
    """
    Remove any environment image setup and go back to black.
    """

    clear_world_background(color=(0.0, 0.0, 0.0, 1.0), strength=1.0)

def clear_collection_objects(collection_name):
    coll = bpy.data.collections.get(collection_name)
    if coll is None:
        return

    for obj in list(coll.objects):
        bpy.data.objects.remove(obj, do_unlink=True)

    bpy.data.collections.remove(coll)


def orthonormal_basis_from_normal(normal):
    n = Vector(normal).normalized()

    # Pick a helper vector not parallel to n
    if abs(n.x) < 0.9:
        helper = Vector((1.0, 0.0, 0.0))
    else:
        helper = Vector((0.0, 1.0, 0.0))

    u = n.cross(helper).normalized()
    v = n.cross(u).normalized()

    return u, v, n


def random_point_on_band(radius, center=(0, 0, 0), band_normal=(0, 0, 1), band_half_angle_deg=12.0):
    """
    Random point on a spherical shell, but concentrated near a plane.
    band_normal sets the pole of the band; the band itself lies around the perpendicular plane.
    """

    u, v, n = orthonormal_basis_from_normal(band_normal)

    az = random.uniform(0.0, 2.0 * math.pi)

    # Elevation above/below the band plane
    max_elev = math.radians(band_half_angle_deg)
    elev = random.uniform(-max_elev, max_elev)

    direction = (
        math.cos(elev) * math.cos(az) * u +
        math.cos(elev) * math.sin(az) * v +
        math.sin(elev) * n
    ).normalized()

    return Vector(center) + radius * direction

def make_emission_material(name="StarEmission", color=(1, 1, 1, 1), strength=5.0):
    """
    Create an emission material for stars.
    """
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)

    mat.use_nodes = True
    nt = mat.node_tree
    nodes = nt.nodes
    links = nt.links

    # Clear existing nodes
    for node in list(nodes):
        nodes.remove(node)

    out = nodes.new(type="ShaderNodeOutputMaterial")
    emission = nodes.new(type="ShaderNodeEmission")

    emission.inputs["Color"].default_value = color
    emission.inputs["Strength"].default_value = strength

    links.new(emission.outputs["Emission"], out.inputs["Surface"])

    return mat

def add_galaxy_band(
    collection_name="GalaxyBand",
    center=(0, 0, 0),
    radius=5000.0,
    num_stars=2500,
    star_radius_min=3.0,
    star_radius_max=10.0,
    emission_strength=20.0,
    band_normal=(0.3, 0.8, 0.5),
    band_half_angle_deg=10.0,
    seed=0,
):
    """
    Add a dense Milky-Way-like band of emissive stars.
    """

    random.seed(seed)

    clear_collection_objects(collection_name)

    coll = bpy.data.collections.new(collection_name)
    bpy.context.scene.collection.children.link(coll)

    star_mats = [
        make_emission_material("GalaxyStarWhite", (1.0, 1.0, 1.0, 1.0), emission_strength),
        make_emission_material("GalaxyStarBlue",  (0.8, 0.9, 1.0, 1.0), emission_strength),
        make_emission_material("GalaxyStarWarm",  (1.0, 0.92, 0.75, 1.0), emission_strength),
    ]

    for i in range(num_stars):
        loc = random_point_on_band(
            radius=radius,
            center=center,
            band_normal=band_normal,
            band_half_angle_deg=band_half_angle_deg,
        )

        r = random.uniform(star_radius_min, star_radius_max)

        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=1,
            radius=r,
            location=loc,
        )

        star = bpy.context.active_object
        star.name = f"GalaxyStar_{i:05d}"

        for old_coll in list(star.users_collection):
            old_coll.objects.unlink(star)
        coll.objects.link(star)

        mat = random.choice(star_mats)
        if star.data.materials:
            star.data.materials[0] = mat
        else:
            star.data.materials.append(mat)

    print(f"Created galaxy band '{collection_name}' with {num_stars} stars.")


def random_point_on_sphere(radius, center=(0, 0, 0)):
    """
    Uniform random point on a sphere.
    """
    u = random.uniform(-1.0, 1.0)
    phi = random.uniform(0.0, 2.0 * math.pi)

    s = math.sqrt(1.0 - u * u)

    x = s * math.cos(phi)
    y = s * math.sin(phi)
    z = u

    return Vector(center) + radius * Vector((x, y, z))


def clear_collection_objects(collection_name):
    """
    Delete a collection and all objects in it.
    """
    coll = bpy.data.collections.get(collection_name)
    if coll is None:
        return

    for obj in list(coll.objects):
        bpy.data.objects.remove(obj, do_unlink=True)

    bpy.data.collections.remove(coll)


def add_star_field(
    collection_name="StarField",
    center=(0, 0, 0),
    radius=5000.0,
    num_stars=300,
    star_radius_min=8.0,
    star_radius_max=20.0,
    emission_strength=4.0,
    seed=0,
):
    """
    Add a shell of emissive stars around the scene.

    center:
        Usually your target point, e.g. (0,0,0)

    radius:
        Distance of the stars from the center.
        Make this much larger than your camera orbit radius.

    num_stars:
        Number of stars to place.

    star_radius_min/max:
        Physical size of the little star spheres.

    emission_strength:
        Brightness of stars.

    seed:
        Random seed for reproducibility.
    """
    random.seed(seed)

    # Remove old star field if present
    clear_star_field(collection_name)

    # Create collection
    coll = bpy.data.collections.new(collection_name)
    bpy.context.scene.collection.children.link(coll)

    # A few star colors
    star_mats = [
        make_emission_material("StarWhite",    (1.0, 1.0, 1.0, 1.0), emission_strength),
        make_emission_material("StarBlue",     (0.85, 0.90, 1.0, 1.0), emission_strength),
        make_emission_material("StarWarm",     (1.0, 0.95, 0.85, 1.0), emission_strength),
    ]

    for i in range(num_stars):
        loc = random_point_on_sphere(radius, center=center)
        r = random.uniform(star_radius_min, star_radius_max)

        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=1,
            radius=r,
            location=loc
        )
        star = bpy.context.active_object
        star.name = f"Star_{i:04d}"

        # Move to star collection
        for old_coll in list(star.users_collection):
            old_coll.objects.unlink(star)
        coll.objects.link(star)

        # Assign a random star material
        mat = random.choice(star_mats)
        if star.data.materials:
            star.data.materials[0] = mat
        else:
            star.data.materials.append(mat)

    print(f"Created {num_stars} stars in collection '{collection_name}'")


def spherical_to_cartesian(radius, polar, azimuth, target=(0, 0, 0), degrees=True):
   """
   Convert spherical coordinates to a Blender world-space position.

   Convention:
     radius  = distance from target
     polar   = angle down from +Z axis
     azimuth = angle around Z axis, measured from +X toward +Y

   So:
     polar = 0     -> camera above target on +Z
     polar = 90    -> camera in XY plane
     azimuth = 0   -> camera on +X side
     azimuth = 90  -> camera on +Y side
   """

   if degrees:
       polar = math.radians(polar)
       azimuth = math.radians(azimuth)

   x = radius * math.sin(polar) * math.cos(azimuth)
   y = radius * math.sin(polar) * math.sin(azimuth)
   z = radius * math.cos(polar)

   return Vector(target) + Vector((x, y, z))


def point_camera_at(camera, target=(0, 0, 0)):
   """
   Rotate camera so it points at target.

   Blender cameras look along their local -Z axis, with local Y as up.
   """

   camera_location = camera.location
   direction = Vector(target) - camera_location

   camera.rotation_euler = direction.to_track_quat('-Z', 'Y').to_euler()


def add_orbit_camera(
   name="OrbitCamera",
   radius=10.0,
   polar=60.0,
   azimuth=0.0,
   target=(0, 0, 0),
   focal_length=50.0,
   clip_start=0.1,
   clip_end=20000.0,
):
   """
   Add a camera at a spherical-coordinate position and point it at target.
   """

   cam_data = bpy.data.cameras.new(name)
   camera = bpy.data.objects.new(name, cam_data)

   bpy.context.collection.objects.link(camera)

   camera.location = spherical_to_cartesian(
       radius=radius,
       polar=polar,
       azimuth=azimuth,
       target=target,
       degrees=True,
   )

   point_camera_at(camera, target)

   cam_data.lens = focal_length
   cam_data.clip_start = clip_start
   cam_data.clip_end = clip_end

   bpy.context.scene.camera = camera

   return camera


def animate_orbit_camera(
   camera,
   frame_start,
   frame_end,
   target=(0, 0, 0),
   radius=10.0,
   polar=60.0,
   azimuth_start=0.0,
   azimuth_end=360.0,
):
   """
   Animate camera along a spherical path.

   This version keeps radius and polar fixed, and sweeps azimuth.
   For example, azimuth 0 -> 360 gives one full orbit.
   """

   num_frames = frame_end - frame_start

   for frame in range(frame_start, frame_end + 1):
       if num_frames == 0:
           t = 0.0
       else:
           t = (frame - frame_start) / num_frames

       azimuth = azimuth_start + t * (azimuth_end - azimuth_start)

       camera.location = spherical_to_cartesian(
           radius=radius,
           polar=polar,
           azimuth=azimuth,
           target=target,
           degrees=True,
       )

       point_camera_at(camera, target)

       camera.keyframe_insert(data_path="location", frame=frame)
       camera.keyframe_insert(data_path="rotation_euler", frame=frame)


def set_camera_interpolation_linear(camera):
   """
   Set camera location and rotation keyframes to LINEAR interpolation.

   Works with older Blender action.fcurves API when available.
   If using a newer Blender version where action.fcurves is gone,
   it silently skips this step.
   """

   if camera.animation_data is None:
       print("Camera has no animation_data; skipping interpolation setup.")
       return

   action = camera.animation_data.action
   if action is None:
       print("Camera has no action; skipping interpolation setup.")
       return

   # Older Blender API
   if hasattr(action, "fcurves"):
       for fcurve in action.fcurves:
           for keyframe in fcurve.keyframe_points:
               keyframe.interpolation = 'LINEAR'
       return

   # Newer Blender API fallback.
   # If fcurves are not directly exposed, just skip.
   print("This Blender version does not expose action.fcurves directly; skipping interpolation setup.")





def animate_parameterized_camera(
   camera,
   frame_start,
   frame_end,
   target=(0, 0, 0),
   r_func=None,
   theta_func=None,
   phi_func=None,
   degrees=True,
):
    """
    Animate a camera using spherical coordinates parameterized in time.

    r_func(t):
       Radius from target.

    theta_func(t):
       Polar angle measured down from +Z.

    phi_func(t):
       Azimuthal angle around Z, measured from +X toward +Y.

    t runs from 0 to 1 over frame_start -> frame_end.
    """

    if r_func is None:
       r_func = lambda t: 10.0

    if theta_func is None:
       theta_func = lambda t: 60.0

    if phi_func is None:
       phi_func = lambda t: 360.0 * t

    num_frames = frame_end - frame_start

    for frame in range(frame_start, frame_end + 1):
        if num_frames == 0:
           t = 0.0
        else:
           t = (frame - frame_start) / num_frames

        r = r_func(t)
        theta = theta_func(t)
        phi = phi_func(t)

        camera.location = spherical_to_cartesian(
           radius=r,
           polar=theta,
           azimuth=phi,
           target=target,
           degrees=degrees,
        )

        point_camera_at(camera, target)

        camera.keyframe_insert(data_path="location", frame=frame)
        camera.keyframe_insert(data_path="rotation_euler", frame=frame)


def setup_image_sequence_render(
   output_dir,
   frame_start,
   frame_end,
   fps=24,
   resolution_x=640,
   resolution_y=360,
   resolution_percentage=100,
   frame_step=1,
   file_prefix="frame_",
   image_format="PNG",
):
    """
    Configure Blender to render animation frames as individual images.

    output_dir:
       Folder where images will be saved.

    file_prefix:
       Output files will look like:
           frame_0001.png
           frame_0002.png
           ...

    frame_step:
       Render every frame if 1.
       Render every 5th frame if 5.
    """

    os.makedirs(output_dir, exist_ok=True)

    scene = bpy.context.scene

    scene.frame_start = frame_start
    scene.frame_end = frame_end
    scene.frame_step = frame_step

    scene.render.fps = fps

    scene.render.resolution_x = resolution_x
    scene.render.resolution_y = resolution_y
    scene.render.resolution_percentage = resolution_percentage

    # Blender automatically appends frame numbers and the image extension.
    scene.render.filepath = os.path.join(output_dir, file_prefix)

    # Faster viewport-like rendering
    # scene.render.engine = 'BLENDER_WORKBENCH'
    # Slower full rendering
    bpy.context.scene.render.engine = 'CYCLES'

    image_settings = scene.render.image_settings

    # Newer Blender versions may distinguish IMAGE vs VIDEO.
    if hasattr(image_settings, "media_type"):
        image_settings.media_type = 'IMAGE'

    image_settings.file_format = image_format
    # Make sure the world/background actually appears in the PNG
    scene.render.film_transparent = False

    if image_format == "PNG":
        image_settings.color_mode = 'RGB'
        image_settings.compression = 15

    print(f"Rendering image sequence to: {output_dir}")
    print(f"Frames: {frame_start} to {frame_end}, step {frame_step}")
    print(f"Resolution: {resolution_x} x {resolution_y} at {resolution_percentage}%")
    print(f"FPS metadata: {fps}")


def render_image_sequence():
   bpy.ops.render.render(animation=True)


#def render_video():
#    """
#    Render the animation.
#    """

#    bpy.ops.render.render(animation=True)
   

def set_render_look(
    background_color=(0.02, 0.02, 0.025),
    exposure=-1.0,
    gamma=1.0,
):
    scene = bpy.context.scene

    scene.render.film_transparent = False

    if scene.world is None:
        scene.world = bpy.data.worlds.new("World")

    scene.world.color = background_color

    scene.view_settings.view_transform = 'Filmic'
    scene.view_settings.look = 'Medium High Contrast'
    scene.view_settings.exposure = exposure
    scene.view_settings.gamma = gamma
   
   
def make_material(
    name="SphereMaterial",
    color=(1.0, 1.0, 1.0, 1.0),
    roughness=0.5,
    metallic=0.0,
):
    """
    Create or update a Principled BSDF material.

    color is RGBA:
        (red, green, blue, alpha)

    Values should usually be between 0 and 1.
    """

    mat = bpy.data.materials.get(name)

    if mat is None:
        mat = bpy.data.materials.new(name)

    mat.use_nodes = True

    bsdf = mat.node_tree.nodes.get("Principled BSDF")

    if bsdf is not None:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = roughness
        bsdf.inputs["Metallic"].default_value = metallic

    mat.diffuse_color = color

    return mat


def set_object_material(obj, material):
    """
    Assign a material to a Blender object.
    """
    if obj.data.materials:
        obj.data.materials[0] = material
    else:
        obj.data.materials.append(material)

def apply_shared_random_oriented_equirectangular_texture_to_spheres(
    image_path,
    material_name="SharedRandomEquirectMaterial",
    sphere_name_patterns=("Sphere", "sphere", "Ball", "ball"),
    roughness=0.8,
    metallic=0.0,
    specular=0.5,
    max_rotation_degrees=360.0,
    do_uv_project=True,
):
    """
    Apply ONE shared material to all sphere-like mesh objects.
    The material uses Object Info -> Random to rotate the texture
    differently for each object.
    """

    mat = make_shared_random_oriented_equirectangular_material(
        image_path=image_path,
        name=material_name,
        roughness=roughness,
        metallic=metallic,
        specular=specular,
        max_rotation_degrees=max_rotation_degrees,
    )

    count = 0

    for obj in bpy.context.scene.objects:
        if obj.type != "MESH":
            continue

        if any(pattern in obj.name for pattern in sphere_name_patterns):
            if do_uv_project:
                sphere_project_uv(obj)

            set_object_material(obj, mat)
            count += 1

    print(f"Applied shared random-oriented texture to {count} sphere-like objects.")

    return mat

def color_spheres(
    color=(0.8, 0.8, 0.8, 1.0),
    material_name="AggregateSphereMaterial",
    sphere_name_patterns=("Sphere", "sphere", "Ball", "ball"),
    roughness=0.5,
    metallic=0.0,
):
    """
    Change the color of all sphere-like mesh objects.

    This searches object names for entries in sphere_name_patterns.
    Change the patterns if your monomers have a different naming scheme.
    """

    mat = make_material(
        name=material_name,
        color=color,
        roughness=roughness,
        metallic=metallic,
    )

    count = 0

    for obj in bpy.context.scene.objects:
        if obj.type != "MESH":
            continue

        if any(pattern in obj.name for pattern in sphere_name_patterns):
            set_object_material(obj, mat)
            count += 1

    print(f"Colored {count} sphere-like objects with material {material_name}")

    return mat


def color_all_meshes(
    color=(0.8, 0.8, 0.8, 1.0),
    material_name="MeshMaterial",
    roughness=0.5,
    metallic=0.0,
):
    """
    Change the color of all mesh objects in the scene.
    Useful if your spheres do not have predictable names.
    """

    mat = make_material(
        name=material_name,
        color=color,
        roughness=roughness,
        metallic=metallic,
    )

    count = 0

    for obj in bpy.context.scene.objects:
        if obj.type == "MESH":
            set_object_material(obj, mat)
            count += 1

    print(f"Colored {count} mesh objects with material {material_name}")

    return mat


def add_point_light(
    name="PointLight",
    location=(0.0, 0.0, 5.0),
    power=500.0,
    color=(1.0, 1.0, 1.0),
    radius=0.0,
):
    """
    Add a point source light.

    power:
        Brightness in Watts.

    radius:
        Larger values make softer shadows.
    """

    light_data = bpy.data.lights.new(name=name, type='POINT')
    light_obj = bpy.data.objects.new(name=name, object_data=light_data)

    bpy.context.collection.objects.link(light_obj)

    light_obj.location = location
    light_data.energy = power
    light_data.color = color
    light_data.shadow_soft_size = radius

    return light_obj

def add_sun_light_from_direction(
    name="SunLight",
    direction=(0.0, 0.0, -1.0),
    power=1.0,
    color=(1.0, 1.0, 1.0),
):
    """
    Add a sun light pointing along a given world-space direction.

    direction is the direction the light rays travel.
    For example:
        direction=(0, 0, -1) means light shines downward from +Z.
        direction=(1, 0, 0) means light shines in the +X direction.
    """

    light_data = bpy.data.lights.new(name=name, type='SUN')
    light_obj = bpy.data.objects.new(name=name, object_data=light_data)

    bpy.context.collection.objects.link(light_obj)

    direction = Vector(direction).normalized()

    # Blender lights emit along their local -Z axis.
    light_obj.rotation_euler = direction.to_track_quat('-Z', 'Y').to_euler()

    light_data.energy = power
    light_data.color = color

    return light_obj


# def add_sun_light(
#     name="SunLight",
#     rotation=(math.radians(45.0), 0.0, math.radians(45.0)),
#     power=2.0,
#     color=(1.0, 1.0, 1.0),
# ):
#     """
#     Add a sun light.

#     Sun lights are directional. Their location does not matter much;
#     the rotation controls the light direction.
#     """

#     light_data = bpy.data.lights.new(name=name, type='SUN')
#     light_obj = bpy.data.objects.new(name=name, object_data=light_data)

#     bpy.context.collection.objects.link(light_obj)

#     light_obj.rotation_euler = rotation
#     light_data.energy = power
#     light_data.color = color

#     return light_obj

def set_dark_background(color=(0.02, 0.02, 0.025, 1.0), strength=1.0):
    """
    Set the rendered world background color in a way that works with Cycles.
    color should be RGBA.
    """

    scene = bpy.context.scene

    if scene.world is None:
        scene.world = bpy.data.worlds.new("World")

    world = scene.world
    world.use_nodes = True

    nt = world.node_tree
    nodes = nt.nodes

    bg = nodes.get("Background")
    if bg is None:
        bg = nodes.new(type="ShaderNodeBackground")

    bg.inputs["Color"].default_value = color
    bg.inputs["Strength"].default_value = strength

    scene.render.film_transparent = False

def clear_lights():
    """
    Delete all light objects in the Blender file.
    """

    for obj in list(bpy.data.objects):
        if obj.type == "LIGHT":
            bpy.data.objects.remove(obj, do_unlink=True)

    # Optional: remove unused light data blocks too
    for light in list(bpy.data.lights):
        if light.users == 0:
            bpy.data.lights.remove(light) 
   
   
def clear_cameras(verbose=True):
    """
    Delete all camera objects in the Blender file.
    """

    camera_objects = [obj for obj in bpy.data.objects if obj.type == "CAMERA"]

    if verbose:
        print(f"Found {len(camera_objects)} camera object(s) to delete.")

    for obj in camera_objects:
        if verbose:
            print(f"Deleting camera object: {obj.name}")
        bpy.data.objects.remove(obj, do_unlink=True)

    unused_cameras = [cam for cam in bpy.data.cameras if cam.users == 0]

    for cam in unused_cameras:
        if verbose:
            print(f"Removing unused camera data: {cam.name}")
        bpy.data.cameras.remove(cam)

    bpy.context.scene.camera = None


def clear_stars():
    clear_collection_objects("StarField")


def clear_galaxy():
    clear_collection_objects("GalaxyBand")

clear_lights()
clear_cameras()
clear_stars()
clear_galaxy()
clear_equirectangular_background()

print("setting background texture")
set_equirectangular_background_visible_only(
    image_path="/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/videos/galaxy.jpg",
    camera_strength=1.0,
    lighting_strength=0.05,
    rotation_z_deg=0.0,
)

# apply_equirectangular_texture_to_spheres(
#     image_path="/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/videos/rock.jpg",
#     material_name="GalaxyDustTexture",
#     roughness=0.9,
#     metallic=0.0,
#     do_uv_project=True,
# )

print("setting sphere texture")
apply_shared_random_oriented_equirectangular_texture_to_spheres(
    image_path="/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/videos/rock.jpg",
    material_name="GalaxyDustShared",
    roughness=0.9,
    metallic=0.0,
    max_rotation_degrees=360.0,
    do_uv_project=False,
)

# Change aggregate color
# color_spheres(
#     color=(0.0, 0.75, 0.0, 1.0),
#     # color=(0.65, 0.65, 0.70, 1.0),
#     material_name="DustAggregateGray",
#     sphere_name_patterns=("Sphere", "sphere", "Ball", "ball"),
#     roughness=0.8,
# )

# Add a soft point light near the camera/front
print("adding point light")
add_point_light(
    name="KeyPointLight",
    location=(0.0, 0.0, 0.0),
    power=800.0,
    color=(1.0, 1, 1),
    radius=0.1,
)


# set_dark_background((0.02, 0.02, 0.025))

sun_direction = Vector((-1.0, -1.0, -1.0)).normalized()

print("adding sun 1")
add_sun_light_from_direction(
    name="SunLightFront",
    direction=sun_direction,
    power=10.0,
)

print("adding sun 2")
add_sun_light_from_direction(
    name="SunLightBack",
    direction=-sun_direction,
    power=1.0,
)


# Suppose your aggregate animation already occupies these Blender frames:
frame_start = 0#bpy.data.scenes[0].frame_start
frame_end = bpy.data.scenes[0].frame_end
print(frame_start,frame_end)

# Point the camera at the aggregate center.
# You can replace this with your aggregate COM if you have it.
target_point = (0.0, 0.0, 0.0)



def r_of_t(t):
   return 1000

def theta_of_t(t):
   return 180*(13.0/51.0)
   # return 60.0 + 15.0 * math.sin(2.0 * math.pi * t)

def phi_of_t(t):
   return 360.0 * t

print('adding orbit camera')
camera = add_orbit_camera(
   name="AggregateOrbitCamera",
   radius=r_of_t(0),
   polar=theta_of_t(0),
   azimuth=phi_of_t(0),
   target=target_point,
   focal_length=70.0,
)

print("animating camera")
animate_parameterized_camera(
   camera=camera,
   frame_start=frame_start,
   frame_end=frame_end,
   target=(0, 0, 0),
   r_func=r_of_t,
   theta_func=theta_of_t,
   phi_func=phi_of_t,
)

print("setting up rendering image sequence")
setup_image_sequence_render(
   # output_dir="/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/videos/aggregate_orbit/",
   output_dir="/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/videos/bigbox_orbit/",
   frame_start=frame_start,
   frame_end=frame_end,
   fps=20,
   resolution_x=1920,
   resolution_y=1080,
   resolution_percentage=100,
   frame_step=1,
   file_prefix="aggregate_",
   image_format="PNG",
)

# set_dark_background((0.0, 0.0, 0.0, 1.0), strength=1.0)
# set_render_look(
#     background_color=(0.0, 0.0, 0.0),
#     exposure=-1.0,
# )

print("rendering image sequence")
render_image_sequence()

# ffmpeg -framerate 20 -pattern_type glob -i "aggregate_*.png" -c:v libx264 -pix_fmt yuv420p bigbox_orbit.mp4