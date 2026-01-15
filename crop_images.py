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
CROP_LEFT   = 500
CROP_BOTTOM = 100
CROP_RIGHT  = 500
CROP_TOP    = 100
# =================================================

os.makedirs(output_dir, exist_ok=True)

valid_exts = {".png", ".jpg", ".jpeg", ".tiff", ".bmp"}

less_crop = ["/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfes-lognormrelax_a-12_N-300_T-3.png","/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/visual_Pfee-lognormrelax_a-12_N-300_T-3.png"]
# images = []
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

        # Compute crop box: (left, upper, right, lower)
        if in_path not in less_crop:
            BOTTOM = CROP_BOTTOM + 100
            TOP = CROP_TOP + 100
        else:
            BOTTOM = CROP_BOTTOM
            TOP = CROP_TOP
        left   = CROP_LEFT
        upper  = TOP
        right  = w - CROP_RIGHT
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



# ---- figure layout ----
fig, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(8, 11),   # tuned for two-column full-width figure
    # constrained_layout=True
)
plt.subplots_adjust(hspace=-0.38, wspace=0.05)

axes = axes.flatten()

for ax, img_path, label in zip(axes, images, labels):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.axis("off")

    # panel label (top-left corner)
    ax.text(
        0.02, 0.95, label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left"
    )

# ---- output ----
# plt.savefig(
#     "porosity_visualization_composite.png",
#     dpi=300,
#     bbox_inches="tight"
# )
plt.savefig(
    f"{folder}porosity_visualization_composite.pdf",
    bbox_inches="tight"
)

plt.close()