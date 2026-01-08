import os
from PIL import Image

# ================= USER SETTINGS =================
input_dir = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals"          # folder containing original images
output_dir = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/data/figures/aggRenders/metricVisuals/cropped" # where cropped images will be saved

# Crop amounts in pixels
# (left, bottom, right, top)
CROP_LEFT   = 500
CROP_BOTTOM = 100
CROP_RIGHT  = 500
CROP_TOP    = 100
# =================================================

os.makedirs(output_dir, exist_ok=True)

valid_exts = {".png", ".jpg", ".jpeg", ".tiff", ".bmp"}

for fname in os.listdir(input_dir):
    name, ext = os.path.splitext(fname)
    ext = ext.lower()

    if ext not in valid_exts:
        continue

    in_path  = os.path.join(input_dir, fname)
    out_path = os.path.join(output_dir, f"{name}_cropped{ext}")

    with Image.open(in_path) as img:
        w, h = img.size

        # Compute crop box: (left, upper, right, lower)
        left   = CROP_LEFT
        upper  = CROP_TOP
        right  = w - CROP_RIGHT
        lower  = h - CROP_BOTTOM

        if left >= right or upper >= lower:
            raise ValueError(
                f"Crop too large for image {fname}: "
                f"resulting size would be invalid."
            )

        cropped = img.crop((left, upper, right, lower))
        cropped.save(out_path)

        print(f"Cropped: {fname} → {out_path}")
