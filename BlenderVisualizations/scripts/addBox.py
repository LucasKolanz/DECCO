import bpy

def add_transparent_box(
    size=(1.0, 1.0, 1.0),
    center=(0.0, 0.0, 0.0),
    color=(0.2, 0.6, 1.0, 0.25),  # RGBA, last value is transparency
    name="Transparent Box"
):
    """
    Add a partially transparent box centered at `center`.

    Parameters
    ----------
    size : tuple
        Box dimensions in x, y, z.
    center : tuple
        Box center location.
    color : tuple
        RGBA color. Alpha < 1 makes it transparent.
    name : str
        Object/material name.
    """

    # Add a unit cube centered at the origin
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=center)
    box = bpy.context.object
    box.name = name

    # Scale cube to requested physical dimensions
    box.dimensions = size

    # Apply scale so dimensions are actually baked into the mesh
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    # Create transparent material
    mat = bpy.data.materials.new(name + " Material")
    mat.diffuse_color = color
    mat.use_nodes = True

    bsdf = mat.node_tree.nodes.get("Principled BSDF")

    if bsdf is not None:
        # Blender 4.x
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = color[3]

        if "Base Color" in bsdf.inputs:
            bsdf.inputs["Base Color"].default_value = color

    # Transparency settings
    mat.blend_method = 'BLEND'
    mat.use_screen_refraction = True
    mat.show_transparent_back = True

    # For Eevee / viewport sorting
    mat.alpha_threshold = 0.01

    box.data.materials.append(mat)

    return box



scaleUp = 1e5

x = 1e-4*scaleUp
y = 1e-4*scaleUp
z = 1e-4*scaleUp

add_transparent_box(
    size=(x,y,z),
    color=(1.0, 0.2, 0.1, 0.2)
)