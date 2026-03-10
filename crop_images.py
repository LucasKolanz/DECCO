import os
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# ================= USER SETTINGS =================
input_dir = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals"          # folder containing original images
output_dir = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/cropped" # where cropped images will be saved
folder = output_dir+'/'
# Crop amounts in pixels
# (left, bottom, right, top)
# CROP_LEFT   = 500
# CROP_BOTTOM = 100
# CROP_RIGHT  = 500
# CROP_TOP    = 100

CROP_LEFT   = 500
CROP_BOTTOM = 100
CROP_RIGHT  = 500
CROP_TOP    = 100

crop0 = [CROP_LEFT,CROP_BOTTOM+50,CROP_RIGHT,CROP_TOP+50]
crop1 = [CROP_LEFT,CROP_BOTTOM,CROP_RIGHT,CROP_TOP]
crop2 = [CROP_LEFT,CROP_BOTTOM+50,CROP_RIGHT,CROP_TOP+75]
# =================================================

os.makedirs(output_dir, exist_ok=True)

valid_exts = {".png", ".jpg", ".jpeg", ".tiff", ".bmp"}

less_crop = ["/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfes-lognormrelax_a-12_N-300_T-3.png","/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfee-lognormrelax_a-12_N-300_T-3.png"]
# images = []

crops = {
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pabc-lognormrelax_a-12_N-300_T-3.png": crop1,
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_PKBM-lognormrelax_a-12_N-300_T-3.png": crop1,
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfee-lognormrelax_a-12_N-300_T-3.png": crop1,
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfes-lognormrelax_a-12_N-300_T-3.png": crop1,
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pch-lognormrelax_a-12_N-300_T-3.png": crop1,
    "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_gcs-lognormrelax_a-12_N-300_T-3.png": crop1,
}

for fname in os.listdir(input_dir):
    name, ext = os.path.splitext(fname)
    ext = ext.lower()

    if ext not in valid_exts:
        continue

    in_path  = os.path.join(input_dir, fname)
    out_path = os.path.join(output_dir, f"{name}_cropped{ext}")
    # images.append(out_path)
    with Image.open(in_path) as img:
        w, h = img.size

        LEFT,BOTTOM,RIGHT,TOP = crops[in_path]
        # BOTTOM = CROP_BOTTOM
        # TOP = CROP_TOP
        left   = LEFT
        upper  = TOP
        right  = w - RIGHT
        lower  = h - BOTTOM

        if left >= right or upper >= lower:
            raise ValueError(
                f"Crop too large for image {fname}: "
                f"resulting size would be invalid."
            )

        cropped = img.crop((left, upper, right, lower))
        cropped.save(out_path)

        print(f"Cropped: {fname} → {out_path}")


images = [
    f"{folder}visual_Pabc-lognormrelax_a-12_N-300_T-3_cropped.png",
    f"{folder}visual_PKBM-lognormrelax_a-12_N-300_T-3_cropped.png",
    f"{folder}visual_Pfee-lognormrelax_a-12_N-300_T-3_cropped.png",
    f"{folder}visual_Pfes-lognormrelax_a-12_N-300_T-3_cropped.png",
    f"{folder}visual_Pch-lognormrelax_a-12_N-300_T-3_cropped.png",
    f"{folder}visual_gcs-lognormrelax_a-12_N-300_T-3_cropped.png",
]

labels = [
    r"(a) $\mathcal{P}_{abc}$",
    r"(b) $\mathcal{P}_{KBM}$",
    r"(c) $\mathcal{P}_{fee}$",
    r"(d) $\mathcal{P}_{fes}$",
    r"(e) $\mathcal{P}_{ch}$",
    r"(f) $\mathcal{P}_{gcs}$",
]

names = [
    r"Equivalent Ellipsoid Porosity",
    r"Gyration Radius Based Porosity",
    r"Fully Enclosing Ellipsoid Porosity",
    r"Fully Enclosing Sphere Porosity",
    r"Convex Hull Porosity",
    r"Geometric Cross Section Porosity",
]

# names = [
#     r"Equivalent Ellipsoid",
#     r"Gyration Radius",
#     r"Enclosing Ellipsoid",
#     r"Enclosing Sphere",
#     r"Convex Hull",
#     r"Geometric Cross Section",
# ]



# ---- figure layout ----
fig, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(8, 11),   # tuned for two-column full-width figure
    # constrained_layout=True
)
# plt.subplots_adjust(hspace=-0.35, wspace=0.0)
plt.subplots_adjust(hspace=-0.00, wspace=-0.1)

axes = axes.flatten()

# heights = [0.05,0.05,0.01,0.01,0.05,0.05]
heights = [0.01,0.01,0.01,0.01,0.01,0.01]
heights = [0.02 for i in heights]

for ax, img_path, label, name, height in zip(axes, images, labels, names, heights):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.7)   # thickness
        spine.set_color("0.1")   # border color
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.axis("off")

    # panel label (top-left corner)
    ax.text(
        0.02, 0.95, label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left"
    )

    # panel label (top-right corner)
    ax.text(
        0.5, height, name,
        transform=ax.transAxes,
        fontsize=11,
        # fontweight="bold",
        va="bottom",
        ha="center"
    )

# ---- output ----
# plt.savefig(
#     "porosity_visualization_composite.png",
#     dpi=300,
#     bbox_inches="tight"
# )
plt.savefig(
    f"{folder}porosity_visualization_composite.pdf",
    bbox_inches="tight",
    dpi=600
)

plt.close()