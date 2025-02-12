import matplotlib.pyplot as plt
import numpy as np
import os
import json
from PIL import Image
from mpl_toolkits.axes_grid1 import ImageGrid

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

plt.rcdefaults()
plt.rcParams['font.size'] = 30

def main():
    fig = plt.figure(figsize=(18, 9))  # Increase height slightly

    # Create ImageGrid with minimal padding
    grid = ImageGrid(fig, 111, nrows_ncols=(1, 4), axes_pad=0.05, share_all=True)

    with open(project_path + "default_files/default_input.json", 'r') as fp:
        input_json = json.load(fp)

    path = input_json["data_directory"]
    image_path = path + "data/figures/aggRenders/"

    image_names = [
        'agg-lognorm_rela_a-0_N-300_T-1000_cropped.png',
        'agg-M__a-0_M-20_N-300_T-1000_cropped.png',
        'agg-M__a-0_M-60_N-300_T-1000_cropped.png',
        'agg-M_a-10_M-100_N-300_T-1000_cropped.png'
    ]

    images = []
    cmaps = []

    for im_name in image_names:
        im = Image.open(image_path + im_name)
        image_array = np.array(im)

        # Determine colormap (grayscale or RGB)
        if image_array.ndim == 2:
            cmaps.append('gray')
        else:
            cmaps.append(None)

        images.append(im)

    # Plot images
    for axe, im, cmap in zip(grid, images, cmaps):
        axe.imshow(im, cmap=cmap, aspect='auto')  # Allow auto aspect ratio
        axe.set_xticks([])  # Hide individual image x-ticks
        axe.set_yticks([])  # Hide individual image y-ticks
        axe.set_frame_on(False)  # Remove border

    # Adjust layout to remove excess whitespace
    fig.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.05)  # Reduce whitespace

    # Create a separate axis for the x-labels
    ax = fig.add_axes([0.1, 0.15, 0.8, 0.05])  # Move x-axis higher
    ax.set_xticks([0.03, 0.35, 0.62, 0.93])  # Adjusted positions for x-axis labels
    ax.set_xticklabels(["M=1", "M=20", "M=60", "M=100"], fontsize=30)

    # Hide y-axis completely
    ax.set_yticks([])  
    ax.yaxis.set_ticklabels([])  
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 5))

    ax.set_xlabel('Fragment Size', fontsize=30)

    # Save with tight bounding box to remove extra whitespace
    plt.savefig(path + 'data/figures/FraggComp_BAPA.png', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
