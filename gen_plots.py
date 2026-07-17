"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file plots the data produced by porosity_FD.py. These are, Porosity abc, Porosity KBM, 
average number of contacts, fractal dimension, bulk density (1-Porosity_KBM), and the final angular momentum. This data is then averaged 
over attempts and saved. The saved data contains the average, the uncertainty, and the number of attempts included in the average.

"""
import os
import sys
import json
import matplotlib.pyplot as plt
import math
import numpy as np

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

import gen_data as gd

styles = ['-','--','-.',':']
colors = ['g','b','orange','r','deeppink','rebeccapurple','y','tab:brown']

def label_from_header(header):

    if header == gd.data_headers[0]:
        return r'$\bm{\mathcal{P}_{abc}}$'
    elif header == gd.data_headers[1]:
        return r'$\bm{\mathcal{P}_{KBM}}$'
    elif header == gd.data_headers[2]:
        return r'$\bm{\mathcal{ANC}}$'
    elif header == gd.data_headers[3]:
        return r'$\bm{\mathcal{D}_{f}}$'
    elif header == gd.data_headers[4]:
        return r'$\bm{\mathcal{A}}$'
    elif header == gd.data_headers[5]:
        return r'$\bm{\mathcal{S}}$'
    elif header == gd.data_headers[6]:
        return r'$\bm{\mathcal{P}_{fee}}$'
    elif header == gd.data_headers[7]:
        return r'$\bm{\mathcal{P}_{fes}}$'
    elif header == gd.data_headers[8]:
        return r'$\bm{\mathcal{P}_{ch}}$'
    elif header == gd.data_headers[9]:
        return r'$\bm{\mathcal{P}_{gcs}}$'
    elif header == gd.data_headers[10]:
        return r'$\bm{r}_{eff}$ ($\bm{\mu m}$)'
        # return r'$\bm{N^{1/3}\left( \frac{S}{\pi \sum_i r_{i}^{2}}} \right)$'
        # return r'$\bm{\langle \sigma \rangle}$'
    else:
        return ""

def plot_with_sigma_bands(ax, x, mean, sigma, label=None, line_kwargs=None, band_labels=True):
    """
    Plot mean and shaded ±1σ (darker) and ±2σ (lighter) bands.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x : array-like
    mean : array-like
    sigma : array-like or scalar
    label : str, label for the mean line
    line_kwargs : dict, extra kwargs for the mean line (linestyle, marker, etc.)
    band_labels : bool, include legend entries for the bands
    """
    x = np.asarray(x)
    mean = np.asarray(mean)
    sigma = np.asarray(sigma) if np.ndim(sigma) else float(sigma)

    # If x might be unsorted (important for log x-scale), sort consistently
    order = np.argsort(x)
    x, mean = x[order], mean[order]
    if np.ndim(sigma):
        sigma = sigma[order]

    lo1, hi1 = mean - sigma, mean + sigma
    lo2, hi2 = mean - 2*sigma, mean + 2*sigma

    # Plot the bands first so they sit behind the mean line
    band2 = ax.fill_between(x, lo2, hi2, alpha=0.20, linewidth=0, label=r"$\pm 2\sigma$" if band_labels else None)
    band1 = ax.fill_between(x, lo1, hi1, alpha=0.35, linewidth=0, label=r"$\pm 1\sigma$" if band_labels else None)

    # Draw the mean on top
    line_kwargs = dict() if line_kwargs is None else dict(line_kwargs)
    mean_line, = ax.plot(x, mean, label=label, zorder=5, **line_kwargs)

    # (Optional) match band facecolor to the line color for a cohesive look
    try:
        base = mean_line.get_color()
        band1.set_facecolor(base)
        band2.set_facecolor(base)
    except Exception:
        pass

    return mean_line, band1, band2


def Tanaka(sizes,initRg,temp):
    Kb = 1.380649e-16; #in erg/K
    rho0 = 2.25 #g/cm^3
    B = 0.4
    Pkbm_c = 0.55

    m = rho0*(4/3)*np.pi*(1e-5)**3
    rhomax = np.sqrt(2)*np.pi*rho0/6
    avgvels = np.sqrt(8*Kb*temp/(np.pi*m))
    Eimp = 0.5*m*avgvels**2
    rhoc = (1-Pkbm_c)*rho0
    Eroll = (1/(4*np.pi))*(3*B*m*avgvels**2*rho0**2)*(rhoc**(-1)-rhomax**(-1))**2##Eimp*1000 # I have no idea what this should be for us
    initV = (4*np.pi/3)*(5/3)**(3/2)*initRg**3
    
    print(f"Eimp: {Eimp}")
    print(f"Eroll: {Eroll}")
    print(f"avgvels: {avgvels}")

    Vols = [initV]
    for i,size in enumerate(sizes[:-1]):
        # print(1/np.sqrt((Vols[i]/(size-1) - m/rhomax)**(-2) +(3*B*Eimp*rho0**2)/(2*np.pi*size*Eroll*m**2)))
        Vfin = size*((1/np.sqrt((Vols[i]/(size-1) - m/rhomax)**(-2) +\
                        (3*B*Eimp*rho0**2)/(2*np.pi*size*Eroll*m**2))) + m/rhomax)

        Vols.append(Vfin)


    return 1-np.array(sizes)*(m/rho0)/np.array(Vols)

def gen_agg_im_plot(show_plots=True,save_plots=False):
    import glob as g
    from PIL import Image
    from mpl_toolkits.axes_grid1 import ImageGrid
    import matplotlib.ticker as plticker
    from matplotlib.patches import FancyArrowPatch

    relative_path = ""
    relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
    project_path = os.path.abspath(relative_path) + '/'

    plt.rcdefaults()
    plt.rcParams['font.size'] = 30
    
    constant = 1.5
    fig,ax = plt.subplots(figsize=(constant*1152*6/500 - .10, constant*1080*3/500 ))

    grid = ImageGrid(fig,111,nrows_ncols=(3,6),axes_pad=0,share_all=True)

    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    # job_group = "constrelax"
    # Title = "Constant"
    job_group = "lognormrelax"
    Title = "Lognormal"

    image_path = path + "data/figures/aggRenders/"
    Nums = [30,100,300]
    # Nums = [30,100]
    temps = [3,10,30,100,300,1000]

    cmaps = []

    images = []
    for N in Nums:
        for t in temps:
            # image = g.glob(image_path+'edited/N{}T{}A*'.format(N,t))[0]
            glob_me = image_path+f'edited/Coloredagg-{job_group}*_a-*_N-{N}_T-{t}_cropped.png'
            # glob_me = image_path+f'/agg-{job_group}*_a-*_N-{N}_T-{t}.png'
            # print(glob_me)
            image = g.glob(glob_me)[0]
            im = Image.open(image)

            # Convert the image to a NumPy array
            image_array = np.array(im)

            # Check if the image is grayscale
            if image_array.ndim == 2:
                cmaps.append('gray')
            else:
                cmaps.append(None)

            images.append(im)

    # print(images)
    # grid[0].imshow(images[0])
    # plt.xticks(np.arange(0,1,step=1/6))
    ax.xaxis.set(ticks=[5,15,25,35,45,55,60],
                ticklabels=['3','10','30','100','300','1000',''],
                # label="x"
                )
    # loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    # ax.xaxis.set_major_locator(loc)
    ax.yaxis.set(ticks=[5,15,25,30],ticklabels=['300','100','30',''])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    index = 0
    for axe,im in zip(grid,images):
        axe.axis('off')
        # axe.imshow(im, cmap=cmaps[index])
        axe.imshow(im, cmap=cmaps[index], zorder=0)

        # index

        # axe.imshow(im,cmap=None)
    # axe.imshow(im, cmap=None)
    ax.set_xlabel('Temp (K)')
    ax.set_ylabel('Number of Particles')


    h = 0.61
    start = (0.35, h)    # figure-fraction coords
    end   = (0.68, h)

    # arrow = FancyArrowPatch(
    #     start, end,
    #     transform=fig.transFigure,
    #     arrowstyle='-|>',
    #     mutation_scale=25,
    #     lw=5, color='black',
    #     zorder=200
    # )
    # fig.add_artist(arrow)

    fig.text(
        0.515, end[1]+0.03,       # slightly above the head
        Title,
        transform=fig.transFigure,
        color='black',
        fontsize=30,
        ha='center'
    )

    if save_plots:
        # plt.savefig(path + f'data/figures/ColoredAggComp_{job_group}.png',dpi=1500)
        plt.savefig(path + f'data/figures/ColoredAggComp_{job_group}.png',dpi=600)
    if show_plots:
        plt.show()

def gen_agg_im_plot_BAPA(show_plots=True,save_plots=False):
    import glob as g
    from PIL import Image
    from mpl_toolkits.axes_grid1 import ImageGrid
    import matplotlib.ticker as plticker
    from matplotlib.patches import FancyArrowPatch

    relative_path = ""
    relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
    project_path = os.path.abspath(relative_path) + '/'

    plt.rcdefaults()
    plt.rcParams['font.size'] = 30
    
    constant = 1.5
    # fig,ax = plt.subplots(figsize=(constant*1080*3/500, constant*1080*3/500 ))
    fig,ax = plt.subplots(figsize=(constant*1152*6/500 - .10, constant*1080*3/500 ))
    fig.delaxes(ax)

    grid = ImageGrid(fig,111,nrows_ncols=(3,3),axes_pad=0,share_all=True)

    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    # job_group = "constrelax"
    # Title = "Constant"
    job_group = "CBAPA"
    Title = "CBAPA"

    image_path = path + "data/figures/aggRenders/"
    # Nums = [300]
    # Nums = [30,100]
    temps = [1000]
    M = [3,5,10,15,20,30,50,60,100]
    C=30

    cmaps = []

    images = []
    # for N in Nums:
    for t in temps:
        for m in M:
            n = m*C
            # image = g.glob(image_path+'edited/N{}T{}A*'.format(N,t))[0]
            glob_me = image_path+f'BAPA/ColoredFragg-{job_group}_a-*_M-{m}_N-{n}_T-{t}.png'
            # glob_me = image_path+f'/agg-{job_group}*_a-*_N-{N}_T-{t}.png'
            # print(glob_me)
            glob_ret = g.glob(glob_me)
            if len(glob_ret) > 0:
                image = g.glob(glob_me)[0]
                im = Image.open(image)

                # Convert the image to a NumPy array
                image_array = np.array(im)

                # Check if the image is grayscale
                if image_array.ndim == 2:
                    cmaps.append('gray')
                else:
                    cmaps.append(None)

                images.append(im)

    # print(images)
    # grid[0].imshow(images[0])
    # plt.xticks(np.arange(0,1,step=1/6))
    # ax.xaxis.set(ticks=[5,15,25]#,35,45,55,60],
    #             # ticklabels=[],
    #             # ticklabels=['3','10','30','100','300','1000',''],
    #             # label="x"
    #             )
    # loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    # ax.xaxis.set_major_locator(loc)
    # ax.yaxis.set(ticks=[5,15,25,30],ticklabels=['300','100','30',''])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    index = 0
    for axe,im in zip(grid,images):
        axe.axis('off')
        # axe.imshow(im, cmap=cmaps[index])
        axe.imshow(im, cmap=cmaps[index], zorder=0)

        # index

        # axe.imshow(im,cmap=None)
    # axe.imshow(im, cmap=None)
    # ax.set_xlabel('Temp (K)')
    # ax.set_ylabel('Number of Particles')


    h = 0.61
    start = (0.35, h)    # figure-fraction coords
    end   = (0.68, h)

    # arrow = FancyArrowPatch(
    #     start, end,
    #     transform=fig.transFigure,
    #     arrowstyle='-|>',
    #     mutation_scale=25,
    #     lw=5, color='black',
    #     zorder=200
    # )
    # fig.add_artist(arrow)

    fig.text(
        0.515, end[1]+0.2,       # slightly above the head
        # 0.515, end[1]+0.03,       # slightly above the head
        Title,
        transform=fig.transFigure,
        color='black',
        fontsize=30,
        ha='center'
    )

    if save_plots:
        # plt.savefig(path + f'data/figures/ColoredAggComp_{job_group}.png',dpi=1500)
        plt.savefig(path + f'data/figures/ColoredFraggComp_{job_group}.png',dpi=600)
    if show_plots:
        plt.show()

def gen_agg_im_plot_paper2(show_plots=True, save_plots=False):
    import os
    import json
    import glob as g
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from PIL import Image

    relative_path = ""
    relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
    project_path = os.path.abspath(relative_path) + '/'

    plt.rcdefaults()
    plt.rcParams['font.size'] = 30
    plt.rcParams.update({
        # 'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    with open(project_path + "default_files/default_input.json", 'r') as fp:
        input_json = json.load(fp)

    path = input_json["data_directory"]
    image_path = path + "data/figures/aggRenders/"

    temps = [1000]
    C = 30
    M_sp = 20

    # -----------------------------
    # Spacing controls
    # -----------------------------
    between_group_hspace = -0.3
    within_group_hspace = -0.2
    within_row_wspace = 0.0

    # -----------------------------
    # Default box settings
    # -----------------------------
    default_box = {
        "x0": 0.02,
        "y0": 0.02,
        "width": 0.97,
        "height": 0.97,
        "label_x": -0.035,
        "label_y": 0.50,
        "linewidth": 1.5,
    }

    simulation_groups = [
        {
            "job_group": "BAPA",
            "title": "\\textbf{CAS}\n$N=300$",
            "X": [3, 5, 10, 15, 20, 30, 50, 60, 75, 100, 150],
            "ncols": 6,
            "box": default_box,
        },
        {
            "job_group": "CBAPA",
            "title": "\\textbf{CNP}\n$C=30$",
            "X": [3, 5, 10, 15, 20, 30, 50, 60, 100],
            "ncols": 4,
            "box": default_box,
        },
        {
            "job_group": "DBAPA",
            "title": "\\textbf{CFS}\n$M=20$",
            "X": [1, 2, 3, 4, 5, 6, 10, 15, 20, 30, 60],
            "ncols": 5,
            "box": default_box,
        },
        {
            "job_group": "BAPAWELD",
            "title": "\\textbf{wCAS}\n$N=300$",
            "X": [3, 15, 100],
            "ncols": 3,
            "box": default_box,
        },
    ]

    def get_box_settings(group_box, default_box):
        box = default_box.copy()

        if group_box is not None:
            box.update(group_box)

        return box

    def get_image_path_and_label(job_group, x, t, C):
        """
        Return the image path and label for one aggregate image.
        """
        if job_group == "CBAPA":
            m = x
            c = C
            n = m * c
            xlab = "M"

        elif job_group == "DBAPA":
            c = x
            m = M_sp
            n = m * c
            xlab = "C"

        else:
            m = x
            n = 300
            xlab = "M"

        glob_me = (
            image_path
            + f"paper2/ColoredFragg-{job_group}_a-*_M-{m}_N-{n}_T-{t}.png"
        )

        glob_ret = g.glob(glob_me)

        if len(glob_ret) == 0:
            print(f"Missing image: {glob_me}")
            return None, None

        return glob_ret[0], f"{xlab}={x}"

    def get_crop_from_white_trim(image_path_i, white_thresh=245, pad=2):
        """
        Return tight crop coords around non-white pixels.

        Coords are PIL-style:
            [x0, y0, x1, y1]

        where x1 and y1 are exclusive.
        """
        img = Image.open(image_path_i).convert("RGB")
        arr = np.array(img)

        nonwhite = np.any(arr < white_thresh, axis=2)

        if not np.any(nonwhite):
            return [0, 0, arr.shape[1], arr.shape[0]]

        ys, xs = np.where(nonwhite)

        x0 = max(0, xs.min() - pad)
        y0 = max(0, ys.min() - pad)
        x1 = min(arr.shape[1], xs.max() + pad + 1)
        y1 = min(arr.shape[0], ys.max() + pad + 1)

        return [x0, y0, x1, y1]

    def get_centered_common_crop(image_paths, white_thresh=245, pad=2):
        """
        Find one crop box for all images.

        Important:
        This crop is centered on the original image center, not centered on
        the union of non-white pixels. This preserves the original centering
        of every aggregate and keeps relative aggregate sizes meaningful.
        """
        if len(image_paths) == 0:
            raise RuntimeError("No image paths were found. Cannot compute crop.")

        max_half_width = 0.0
        max_half_height = 0.0

        image_width = None
        image_height = None

        for image_path_i in image_paths:
            img = Image.open(image_path_i).convert("RGB")
            arr = np.array(img)

            h, w = arr.shape[:2]

            if image_width is None:
                image_width = w
                image_height = h
            else:
                if w != image_width or h != image_height:
                    raise ValueError(
                        "Images are not all the same original size. "
                        "A single common crop is only safe if they all have "
                        "the same width and height."
                    )

            cx = w / 2.0
            cy = h / 2.0

            crop_coords = get_crop_from_white_trim(
                image_path_i,
                white_thresh=white_thresh,
                pad=pad
            )

            x0, y0, x1, y1 = crop_coords

            max_half_width = max(
                max_half_width,
                cx - x0,
                x1 - cx
            )

            max_half_height = max(
                max_half_height,
                cy - y0,
                y1 - cy
            )

        x0_common = int(np.floor(image_width / 2.0 - max_half_width))
        x1_common = int(np.ceil(image_width / 2.0 + max_half_width))
        y0_common = int(np.floor(image_height / 2.0 - max_half_height))
        y1_common = int(np.ceil(image_height / 2.0 + max_half_height))

        # Keep crop inside the original image.
        x0_common = max(0, x0_common)
        y0_common = max(0, y0_common)
        x1_common = min(image_width, x1_common)
        y1_common = min(image_height, y1_common)

        return [x0_common, y0_common, x1_common, y1_common]

    def load_and_crop(image_path_i, crop_coords):
        """
        Load one image and crop it using common crop coords.
        """
        img = Image.open(image_path_i).convert("RGB")

        crop_coords = [int(v) for v in crop_coords]
        return img.crop(tuple(crop_coords))

    def load_images_from_paths(image_full_paths, labels, crop_coords):
        images = []
        cmaps = []

        for image_path_i in image_full_paths:
            im = load_and_crop(image_path_i, crop_coords)
            image_array = np.array(im)

            if image_array.ndim == 2:
                cmaps.append("gray")
            else:
                cmaps.append(None)

            images.append(im)

        return images, cmaps, labels

    # -----------------------------
    # First collect image paths and labels
    # -----------------------------
    all_image_paths = []
    group_path_data = []

    for group in simulation_groups:
        image_full_paths = []
        labels = []

        for t in temps:
            for x in group["X"]:
                image, label = get_image_path_and_label(
                    group["job_group"],
                    x,
                    t,
                    C
                )

                if image is None:
                    continue

                image_full_paths.append(image)
                labels.append(label)
                all_image_paths.append(image)

        group_path_data.append({
            "group": group,
            "image_full_paths": image_full_paths,
            "labels": labels,
        })

    # -----------------------------
    # One centered common crop for all images
    # -----------------------------
    common_crop_coords = get_centered_common_crop(all_image_paths)

    print("Common crop coords:", common_crop_coords)

    # -----------------------------
    # Load all images and compute layout
    # -----------------------------
    group_data = []
    max_cols = 0

    for item in group_path_data:
        group = item["group"]
        image_full_paths = item["image_full_paths"]
        labels = item["labels"]

        images, cmaps, labels = load_images_from_paths(
            image_full_paths,
            labels,
            common_crop_coords
        )

        nimgs = len(images)
        ncols = group["ncols"]

        if nimgs == 0:
            print(f"No images found for group {group['title']}")
            nrows_group = 1
        else:
            nrows_group = math.ceil(nimgs / ncols)

        box = get_box_settings(group.get("box", None), default_box)

        group_data.append({
            "title": group["title"],
            "images": images,
            "cmaps": cmaps,
            "labels": labels,
            "ncols": ncols,
            "nrows_group": nrows_group,
            "box": box,
        })

        max_cols = max(max_cols, ncols)

    # -----------------------------
    # Figure layout
    # -----------------------------
    constant = 1.5
    total_group_rows = sum(data["nrows_group"] for data in group_data)

    fig_width = constant * 1.6 * max_cols
    fig_height = constant * 1.6 * total_group_rows

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Explicit spacer rows make the vertical gaps equal even when
    # group heights are different.
    group_height_ratios = [data["nrows_group"] for data in group_data]

    spacer_height = 0.04  # tune this instead of between_group_hspace
    

    outer_height_ratios = []
    for i, ratio in enumerate(group_height_ratios):
        outer_height_ratios.append(ratio)

        # Add a spacer row after every group except the last.
        if i < len(group_height_ratios) - 1:
            outer_height_ratios.append(spacer_height)

    outer_grid = fig.add_gridspec(
        nrows=len(outer_height_ratios),
        ncols=1,
        height_ratios=outer_height_ratios,
        hspace=0.0
    )

    # Make room on the left for group names.
    fig.subplots_adjust(
        left=0.08,
        right=0.995,
        top=0.98,
        bottom=0.02
    )
    # fig.subplots_adjust(
    #     left=0.14,
    #     right=0.98,
    #     top=0.98,
    #     bottom=0.02
    # )

    axes = []
    group_box_axes = []

    # -----------------------------
    # Plot each group
    # -----------------------------
    label_y_offset = [[-0.1,-0.1],[-0.1,0.05],[-0.1,-0.1],[-0.1]]
    group_y_shift = [-0.01,0.005,-0.01,-0.01] 
    for group_i, data in enumerate(group_data):
        grid_row = 2 * group_i
        images = data["images"]
        cmaps = data["cmaps"]
        labels = data["labels"]
        ncols = data["ncols"]
        nrows_group = data["nrows_group"]
        title = data["title"]
        box = data["box"]

        group_axes = []


        # -------------------------------------------------
        # Background axes for this whole group.
        # The rectangle and group label are drawn on this.
        # -------------------------------------------------
        box_ax = fig.add_subplot(outer_grid[grid_row, 0])
        box_ax.set_zorder(100)
        box_ax.patch.set_alpha(0.0)
        box_ax.set_xticks([])
        box_ax.set_yticks([])
        box_ax.set_frame_on(False)


        box_ax.add_patch(
            Rectangle(
                (box["x0"], box["y0"]),
                box["width"],
                box["height"],
                transform=box_ax.transAxes,
                fill=False,
                edgecolor="black",
                linewidth=box["linewidth"],
                clip_on=False,
                zorder=1
            )
        )

        box_ax.text(
            box["label_x"],
            box["label_y"],
            title,
            transform=box_ax.transAxes,
            ha="center",
            va="center",
            fontsize=24,
            rotation=90,
            clip_on=False
        )

        group_box_axes.append(box_ax)

        # -------------------------------------------------
        # Subgrid inside this group.
        # This controls row spacing within the group.
        # -------------------------------------------------
        # First make a padded image area inside the group box.
        # Increasing image_top_pad shifts the image rows downward. 
        image_area_grid = outer_grid[grid_row, 0].subgridspec(
            nrows=3,
            ncols=1,
            height_ratios=[
                0,
                1.0,
                0
            ],
            hspace=0.0
        )

        # Then put the actual image rows only in the middle area.
        group_grid = image_area_grid[1, 0].subgridspec(
            nrows=nrows_group,
            ncols=1,
            hspace=within_group_hspace
        )

        for local_row in range(nrows_group):
            start_idx = local_row * ncols
            end_idx = min(start_idx + ncols, len(images))

            row_images = images[start_idx:end_idx]
            row_cmaps = cmaps[start_idx:end_idx]
            row_labels = labels[start_idx:end_idx]

            n_this_row = len(row_images)

            if n_this_row == 0:
                continue

            # Each row has the same total possible width.
            # Incomplete rows are centered.
            inner_grid = group_grid[local_row, 0].subgridspec(
                nrows=1,
                ncols=2 * max_cols,
                wspace=within_row_wspace
            )

            # Center row in half-image increments.
            start_col = max_cols - n_this_row

            for j, im in enumerate(row_images):
                col0 = start_col + 2 * j
                col1 = col0 + 2

                ax = fig.add_subplot(inner_grid[0, col0:col1])
                image_y_shift = group_y_shift[group_i]

                pos = ax.get_position()
                ax.set_position([
                    pos.x0,
                    pos.y0 + image_y_shift,
                    pos.width,
                    pos.height
                ])
                ax.set_zorder(2)

                ax.imshow(im, cmap=row_cmaps[j])
                ax.set_aspect("equal", adjustable="box")
                ax.set_anchor("C")
                ax.axis("off")

                ax.text(
                    0.5,
                    0.98+label_y_offset[group_i][local_row],
                    row_labels[j],
                    transform=ax.transAxes,
                    ha="center",
                    va="top",
                    fontsize=14
                )

                group_axes.append(ax)

        axes.append(group_axes)

    if save_plots:
        plt.savefig(
            path + "data/figures/ColoredFraggComp_all_groups.png",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.05
        )

    if show_plots:
        plt.show()

    return fig, axes

# def gen_agg_im_plot_paper2(show_plots=True, save_plots=False):
#     import os
#     import json
#     import glob as g
#     import math
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from matplotlib.patches import Rectangle
#     from PIL import Image

#     relative_path = ""
#     relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
#     project_path = os.path.abspath(relative_path) + '/'

#     plt.rcdefaults()
#     plt.rcParams['font.size'] = 30

#     with open(project_path + "default_files/default_input.json", 'r') as fp:
#         input_json = json.load(fp)

#     path = input_json["data_directory"]
#     image_path = path + "data/figures/aggRenders/"

#     temps = [1000]
#     C = 30
#     M_sp = 20

#     # -----------------------------
#     # Spacing controls
#     # -----------------------------
#     between_group_hspace = 0.22   # spacing between the group boxes
#     within_group_hspace = 0.0    # spacing between rows inside each group
#     within_row_wspace = 0.0      # spacing between images in the same row
#     # within_row_wspace = 0.001      # spacing between images in the same row

#     # -----------------------------
#     # Default box settings
#     # Used if a group does not specify one of these values
#     # -----------------------------
#     default_box = {
#         "x0": 0.02,
#         "y0": -0.03,
#         "width": 0.97,
#         "height": 1.2,
#         "label_x": -0.035,
#         "label_y": 0.50,
#         "linewidth": 1.5,
#     }

#     simulation_groups = [
#         {
#             "job_group": "BAPA",
#             "title": "CAS",
#             "X": [3, 5, 10, 15, 20, 30, 50, 60, 75, 100, 150],
#             "ncols": 6,
#             "box": default_box,
#             # "box": {
#             #     "x0": 0.02,
#             #     "y0": 0.04,
#             #     "width": 0.96,
#             #     "height": 0.98,
#             #     "label_x": -0.035,
#             #     "label_y": 0.50,
#             #     "linewidth": 1.5,
#             # },
#         },
#         {
#             "job_group": "CBAPA",
#             "title": "CNP",
#             "X": [3, 5, 10, 15, 20, 30, 50, 60, 100],
#             "ncols": 5,
#             "box": default_box,
#             # "box": {
#             #     "x0": 0.02,
#             #     "y0": 0.04,
#             #     "width": 0.96,
#             #     "height": 0.98,
#             #     "label_x": -0.035,
#             #     "label_y": 0.50,
#             #     "linewidth": 1.5,
#             # },
#         },
#         {
#             "job_group": "DBAPA",
#             "title": "CFS",
#             "M": [20],
#             "X": [1,2,3,4,5,6,10,15,20,30,60],
#             "ncols": 6,
#             "box": default_box,
#             # "box": {
#             #     "x0": 0.02,
#             #     "y0": 0.04,
#             #     "width": 0.96,
#             #     "height": 0.98,
#             #     "label_x": -0.035,
#             #     "label_y": 0.50,
#             #     "linewidth": 1.5,
#             # },
#         },
#         {
#             "job_group": "BAPAWELD",
#             "title": "wCAS",
#             "X": [3, 15, 100],
#             "ncols": 3,
#             "box": default_box,
#             # "box": {
#             #     "x0": 0.02,
#             #     "y0": 0.04,
#             #     "width": 0.96,
#             #     "height": 0.98,
#             #     "label_x": -0.035,
#             #     "label_y": 0.50,
#             #     "linewidth": 1.5,
#             # },
#         },
#     ]

#     def get_box_settings(group_box, default_box):
#         """
#         Merge user-specified box settings with defaults.
#         This lets you omit values from individual groups if you want.
#         """
#         box = default_box.copy()

#         if group_box is not None:
#             box.update(group_box)

#         return box

#     def get_group_crop(job_group, X_values,temps,C):
#         # images = []
#         # cmaps = []
#         # labels = []
#         min_crop_coords = [np.inf,np.inf,-np.inf,-np.inf]
#         # image_full_paths = []

#         for t in temps:
#             for x in X_values:
#                 if job_group == "CBAPA":
#                     m = x
#                     c = C
#                     n = m * c
#                     # xlab = "M"
#                 elif job_group == "DBAPA":
#                     c = x
#                     m = M_sp
#                     n = m * c
#                     # xlab = "C"
#                 else:
#                     m = x
#                     n = 300
#                     # xlab = "M"

#                 glob_me = (
#                     image_path
#                     + f'paper2/ColoredFragg-{job_group}_a-*_M-{m}_N-{n}_T-{t}.png'
#                 )

#                 glob_ret = g.glob(glob_me)

#                 if len(glob_ret) == 0:
#                     print(f"Missing image: {glob_me}")
#                     continue

#                 image = glob_ret[0]
#                 crop_coords = get_crop_from_white_trim(image)
#                 # print(crop_coords)
#                 # print(min_crop_coords)
#                 min_crop_coords = [
#                         np.minimum([crop_coords[0]],[min_crop_coords[0]])[0],
#                         np.minimum([crop_coords[1]],[min_crop_coords[1]])[0],
#                         np.maximum([crop_coords[2]],[min_crop_coords[2]])[0],
#                         np.maximum([crop_coords[3]],[min_crop_coords[3]])[0]]

#         return min_crop_coords

#     def load_group_images(job_group, X_values, temps, C, crop_coords):
#         images = []
#         cmaps = []
#         labels = []
#         # min_crop_coords = [np.inf,np.inf,-np.inf,-np.inf]
#         image_full_paths = []

#         for t in temps:
#             for x in X_values:
#                 if job_group == "CBAPA":
#                     m = x
#                     c = C
#                     n = m * c
#                     xlab = "M"
#                 elif job_group == "DBAPA":
#                     c = x
#                     m = M_sp
#                     n = m * c
#                     xlab = "C"
#                 else:
#                     m = x
#                     n = 300
#                     xlab = "M"

#                 glob_me = (
#                     image_path
#                     + f'paper2/ColoredFragg-{job_group}_a-*_M-{m}_N-{n}_T-{t}.png'
#                 )

#                 glob_ret = g.glob(glob_me)

#                 if len(glob_ret) == 0:
#                     print(f"Missing image: {glob_me}")
#                     continue

#                 image = glob_ret[0]
#                 # crop_coords = get_crop_from_white_trim(image)
#                 # # print(crop_coords)
#                 # # print(min_crop_coords)
#                 # min_crop_coords = [
#                 #         np.minimum([crop_coords[0]],[min_crop_coords[0]])[0],
#                 #         np.minimum([crop_coords[1]],[min_crop_coords[1]])[0],
#                 #         np.maximum([crop_coords[2]],[min_crop_coords[2]])[0],
#                 #         np.maximum([crop_coords[3]],[min_crop_coords[3]])[0]]

#                 image_full_paths.append(image)

#                 labels.append(f"{xlab}={x}")
        
#         for image in image_full_paths:
#             im = load_and_crop(image, crop_coords)
#             image_array = np.array(im)

#             if image_array.ndim == 2:
#                 cmaps.append('gray')
#             else:
#                 cmaps.append(None)
#             images.append(im)

#         return images, cmaps, labels

#     # -----------------------------
#     # Load all images and compute layout
#     # -----------------------------
#     group_data = []
#     max_cols = 0

#     min_crop_coords = [np.inf,np.inf,-np.inf,-np.inf]
#     for group in simulation_groups:
#         crop_coords = get_group_crop(group["job_group"],group["X"],temps,C)
#         min_crop_coords = [
#             np.minimum([crop_coords[0]],[min_crop_coords[0]])[0],
#             np.minimum([crop_coords[1]],[min_crop_coords[1]])[0],
#             np.maximum([crop_coords[2]],[min_crop_coords[2]])[0],
#             np.maximum([crop_coords[3]],[min_crop_coords[3]])[0]]

#     for group in simulation_groups:
#         images, cmaps, labels = load_group_images(
#             group["job_group"],
#             group["X"],
#             temps,
#             C,
#             min_crop_coords
#         )

#         nimgs = len(images)
#         ncols = group["ncols"]

#         if nimgs == 0:
#             print(f"No images found for group {group['title']}")
#             nrows_group = 1
#         else:
#             nrows_group = math.ceil(nimgs / ncols)

#         box = get_box_settings(group.get("box", None), default_box)

#         group_data.append({
#             "title": group["title"],
#             "images": images,
#             "cmaps": cmaps,
#             "labels": labels,
#             "ncols": ncols,
#             "nrows_group": nrows_group,
#             "box": box,
#         })

#         max_cols = max(max_cols, ncols)

#     # -----------------------------
#     # Figure layout
#     # -----------------------------
#     constant = 1.5
#     total_group_rows = sum(data["nrows_group"] for data in group_data)

#     fig_width = constant * 1.6 * max_cols
#     fig_height = constant * 1.6 * total_group_rows

#     fig = plt.figure(figsize=(fig_width, fig_height))

#     # One outer row per group.
#     # This lets between_group_hspace only affect spacing between groups.
#     outer_grid = fig.add_gridspec(
#         nrows=len(group_data),
#         ncols=1,
#         hspace=between_group_hspace
#     )

#     # Make room on the left for group names.
#     fig.subplots_adjust(
#         left=0.14,
#         right=0.98,
#         top=0.98,
#         bottom=0.02
#     )

#     axes = []
#     group_box_axes = []

#     # -----------------------------
#     # Plot each group
#     # -----------------------------
#     for group_i, data in enumerate(group_data):
#         images = data["images"]
#         cmaps = data["cmaps"]
#         labels = data["labels"]
#         ncols = data["ncols"]
#         nrows_group = data["nrows_group"]
#         title = data["title"]
#         box = data["box"]

#         group_axes = []

#         # -------------------------------------------------
#         # Background axes for this whole group.
#         # The rectangle and group label are drawn on this.
#         # -------------------------------------------------
#         box_ax = fig.add_subplot(outer_grid[group_i, 0])
#         box_ax.set_zorder(100)
#         box_ax.patch.set_alpha(0.0)
#         box_ax.set_xticks([])
#         box_ax.set_yticks([])
#         box_ax.set_frame_on(False)

#         box_ax.add_patch(
#             Rectangle(
#                 (box["x0"], box["y0"]),
#                 box["width"],
#                 box["height"],
#                 transform=box_ax.transAxes,
#                 fill=False,
#                 edgecolor="black",
#                 linewidth=box["linewidth"],
#                 clip_on=False,
#                 zorder=1
#             )
#         )

#         box_ax.text(
#             box["label_x"],
#             box["label_y"],
#             title,
#             transform=box_ax.transAxes,
#             ha="right",
#             va="center",
#             fontsize=24,
#             rotation=90,
#             clip_on=False
#         )

#         group_box_axes.append(box_ax)

#         # -------------------------------------------------
#         # Subgrid inside this group.
#         # This controls row spacing within the group.
#         # -------------------------------------------------
#         group_grid = outer_grid[group_i, 0].subgridspec(
#             nrows=nrows_group,
#             ncols=1,
#             hspace=within_group_hspace
#         )

#         for local_row in range(nrows_group):
#             start_idx = local_row * ncols
#             end_idx = min(start_idx + ncols, len(images))

#             row_images = images[start_idx:end_idx]
#             row_cmaps = cmaps[start_idx:end_idx]
#             row_labels = labels[start_idx:end_idx]

#             n_this_row = len(row_images)

#             # Skip empty rows, which can happen if no images were found.
#             if n_this_row == 0:
#                 continue

#             # This makes every row use the same possible total width.
#             # Incomplete rows are centered.
#             inner_grid = group_grid[local_row, 0].subgridspec(
#                 nrows=1,
#                 ncols=2 * max_cols,
#                 wspace=within_row_wspace
#             )

#             # Center row in half-image increments.
#             start_col = max_cols - n_this_row

#             for j, im in enumerate(row_images):
#                 col0 = start_col + 2 * j
#                 col1 = col0 + 2

#                 ax = fig.add_subplot(inner_grid[0, col0:col1])
#                 ax.set_zorder(2)

#                 ax.imshow(im, cmap=row_cmaps[j])
#                 ax.axis("off")
#                 ax.set_title(row_labels[j], fontsize=14, pad=0)

#                 group_axes.append(ax)

#         axes.append(group_axes)

#     if save_plots:
#         plt.savefig(
#             path + 'data/figures/ColoredFraggComp_all_groups.png',
#             dpi=600,
#             bbox_inches="tight",
#             pad_inches=0.05
#         )

#     if show_plots:
#         plt.show()

#     return fig, axes

def gen_relax_vs_tense_seqstick_plots(distribution,show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    tense_data_prefolder = path + 'jobs/SeqStickLognorm_'
    relax_data_prefolder = path + 'jobs/SeqStickLognormrelax_'

    dataset_name = tense_data_prefolder.split("/")[-1]

    figure_folder = path+'data/figures/'


    temps = [1000]
    # temps = [3,10]
    Nums = [300]

    
    
    attempts = [i for i in range(30)]


    # requested_data_headers = gd.data_headers[:2]
    requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]



    tense_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    relax_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for n_i,n in enumerate(Nums):
            for t_i,t in enumerate(temps):
                for data_prefolder in [tense_data_prefolder,relax_data_prefolder]:
                    rel = ""
                    if data_prefolder == tense_data_prefolder:
                        raw_data = tense_raw_data
                    elif data_prefolder == relax_data_prefolder:
                        raw_data = relax_raw_data
                        rel="relax_"

                    folder = f"{data_prefolder}{a}/N_{n}/"
                    print(folder)

                    if os.path.exists(folder+f"{rel}job_data.csv"):
                        with open(folder+f"{rel}job_data.csv",'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if n not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(n)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        print(existing_headers_for_size)
                        print(existing_values_for_size)

                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                print(raw_data.size)
                                print(existing_headers_for_size.index(header))
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)

    tense_avg_data = np.nanmean(tense_raw_data,axis=1)
    tense_std_data = np.nanstd(tense_raw_data,axis=1)
    tense_num_data = np.count_nonzero(~np.isnan(tense_raw_data),axis=1)
    tense_err_data = tense_std_data/np.sqrt(tense_num_data)

    relax_avg_data = np.nanmean(relax_raw_data,axis=1)
    relax_std_data = np.nanstd(relax_raw_data,axis=1)
    relax_num_data = np.count_nonzero(~np.isnan(relax_raw_data),axis=1)
    relax_err_data = relax_std_data/np.sqrt(relax_num_data)


    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Tense data has {} nan values".format(np.count_nonzero(np.isnan(tense_avg_data))))
    print("Relax data has {} nan values".format(np.count_nonzero(np.isnan(relax_avg_data))))
    
    print("For tense data:")
    for h_i,header in enumerate(requested_data_headers):
        print(f"\t{distribution} {header}: {tense_avg_data[h_i]} +- {tense_err_data[h_i]}")

    print("For relax data:")
    for h_i,header in enumerate(requested_data_headers):
        print(f"\t{distribution} {header}: {relax_avg_data[h_i]} +- {relax_err_data[h_i]}")


def gen_relax_vs_tense_BPCA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    tense_data_prefolder = path + 'jobsCosine/lognorm_'
    relax_data_prefolder = path + 'jobsCosine/lognormrelax_'

    dataset_name = tense_data_prefolder.split("/")[-1]

    figure_folder = path+'data/figures/'


    temps = [3,10,30,100,300,1000]
    # temps = [3,10]
    Nums = [30,100,300]

    
    
    attempts = [i for i in range(30)]


    # requested_data_headers = gd.data_headers[:2]
    requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]



    tense_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    relax_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for n_i,n in enumerate(Nums):
            for t_i,t in enumerate(temps):
                for data_prefolder in [tense_data_prefolder,relax_data_prefolder]:
                    rel = ""
                    if data_prefolder == tense_data_prefolder:
                        raw_data = tense_raw_data
                    elif data_prefolder == relax_data_prefolder:
                        raw_data = relax_raw_data
                        rel="relax_"

                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    print(folder)

                    if os.path.exists(folder+f"{rel}job_data.csv"):
                        with open(folder+f"{rel}job_data.csv",'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if n not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(n)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(gd.data_headers):
                            if header in existing_headers_for_size:

                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)

    tense_avg_data = np.nanmean(tense_raw_data,axis=1)
    tense_std_data = np.nanstd(tense_raw_data,axis=1)
    tense_num_data = np.count_nonzero(~np.isnan(tense_raw_data),axis=1)
    tense_err_data = tense_std_data/np.sqrt(tense_num_data)

    relax_avg_data = np.nanmean(relax_raw_data,axis=1)
    relax_std_data = np.nanstd(relax_raw_data,axis=1)
    relax_num_data = np.count_nonzero(~np.isnan(relax_raw_data),axis=1)
    relax_err_data = relax_std_data/np.sqrt(relax_num_data)


    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Tense data has {} nan values".format(np.count_nonzero(np.isnan(tense_avg_data))))
    print("RElax data has {} nan values".format(np.count_nonzero(np.isnan(relax_avg_data))))
    

    
    length = len(temps)


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):
        for n_i,n in enumerate(Nums):


            fig,ax = plt.subplots()

            print(tense_avg_data[h_i,n_i,:])
            print(relax_avg_data[h_i,n_i,:])

            ax.errorbar(temps,tense_avg_data[h_i,n_i,:],yerr=tense_err_data[h_i,n_i,:],\
                    label=f"tense N={n}",\
                    linestyle=styles[h_i],marker='.',markersize=10,zorder=5)
            ax.errorbar(temps,relax_avg_data[h_i,n_i,:],yerr=relax_err_data[h_i,n_i,:],\
                    label=f"relax N={n}",\
                    linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

            # if include_totals:
            #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
            #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            ax.set_xlabel('Temp K')
            ax.set_ylabel(header)
            # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
            ax.set_xscale('log')
            # if i == 1:
            fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
            plt.tight_layout()
            if save_plots:
                plt.savefig("{}{}_{}_tenseVsRelax.png".format(figure_folder,dataset_name,header))
            if show_plots:
                plt.show() 

def gen_Asym_BAPA_numbers():
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolder = path + 'jobs/AsymBAPA_'

    dataset_name = data_prefolder.split("/")[-1]


    temps = [1000]
    # temps = [3,10]]
    # M = [1,3,5,10,15,20,30,50,60,100]
    N = 300
    
    
    attempts = [i for i in range(30)]

    requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]
    


    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for t_i,t in enumerate(temps):
            folder = f"{data_prefolder}{a}/N_{N}/T_{t}/"
            if os.path.exists(folder+"job_data.csv"):
                with open(folder+"job_data.csv",'r') as fp:
                    existing_data = fp.readlines()

                existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                #even though the data can have other sizes in it, 
                #we only want the data of size n
                
                n = u.find_max_index(folder)

                if n not in existing_sizes:
                    print(f"ERROR: Data of size {n} does not exist for {folder}.")
                    continue
                index = existing_sizes.index(n)*4
                existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                
                for h_i,header in enumerate(requested_data_headers):
                    if header in existing_headers_for_size:
                        raw_data[h_i,a_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)

    avg_data = np.nanmean(raw_data,axis=1)
    std_data = np.nanstd(raw_data,axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data = std_data/np.sqrt(num_data)

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(avg_data))))
    


    length = len(temps)


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):
        for t_i,t in enumerate(temps):

            print(f"{header} = {avg_data[h_i,t_i]} +- {err_data[h_i,t_i]}")
    print(f"Average over {num_data}.")

def gen_BAPA_eff_rad(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    relax = False

    data_prefolders = [path + 'jobs/BAPA_', path + 'jobs/CBAPA_', path + 'jobs/BAPAWELD_']


    data_prefolder = data_prefolders[0]
    dataset_name = data_prefolder.split("/")[-1]
    figure_folder = path+'data/figures/'

    temps = [1000]
    Nums = [300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150]
    
    attempts = [i for i in range(30)]


    # requested_data_headers = [1,1,0,1,1,0,0,0,0,1,1]
    requested_data_headers = [gd.data_headers[10]]
    print(requested_data_headers)
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            n = 300
            size = n
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)


    avg_data_BAPA = np.nanmean(raw_data,axis=1)
    std_data_BAPA = np.nanstd(raw_data,axis=1)
    num_data_BAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPA = std_data_BAPA/np.sqrt(num_data_BAPA)
    


    data_prefolder = data_prefolders[1]
    dataset_name = data_prefolder.split("/")[-1]

    temps = [1000]
    C = 30
    # M = [1,3,5,10,15,20,30,50,60,100,150]
    attempts = [i for i in range(30)]

    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            n = C*m
            size = n
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)

    avg_data_CBAPA = np.nanmean(raw_data,axis=1)
    std_data_CBAPA = np.nanstd(raw_data,axis=1)
    num_data_CBAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_CBAPA = std_data_CBAPA/np.sqrt(num_data_CBAPA)


    data_prefolder = data_prefolders[2]
    dataset_name = data_prefolder.split("/")[-1]

    temps = [1000]
    # M = [3,100]
    n = 300
    attempts = [i for i in range(30)]

    # requested_data_headers = [1,1,0,1,1,0,0,0,0,1]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
                            # raw_data[h_i,a_i,m_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

    avg_data_BAPAWELD = np.nanmean(raw_data,axis=1)
    std_data_BAPAWELD = np.nanstd(raw_data,axis=1)
    num_data_BAPAWELD = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPAWELD = std_data_BAPAWELD/np.sqrt(num_data_BAPAWELD)

    for m_i,m in enumerate(M):
        if m in [3,15,100]:
            for h_i,header in enumerate(requested_data_headers):
                # print(f"BAPA M={m} {header}    : {avg_data_BAPA[h_i,m_i,0]}+-{err_data_BAPA[h_i,m_i,0]} ({num_data_BAPA[h_i,m_i,0]})")
                # print(f"BAPAWELD M={m} {header}: {avg_data_BAPAWELD[h_i,m_i,0]}+-{err_data_BAPAWELD[h_i,m_i,0]} ({num_data_BAPAWELD[h_i,m_i,0]})")
                min_val = np.min([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])
                min_unc = [err_data_BAPA[h_i,m_i,0],err_data_BAPAWELD[h_i,m_i,0]][np.argmin([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])]
                max_val = np.max([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])
                max_unc = [err_data_BAPA[h_i,m_i,0],err_data_BAPAWELD[h_i,m_i,0]][np.argmax([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])]
                # if min_val + min_unc >= max_val-max_unc:
                #     print(f"{header} agrees")
                # else:
                #     print(f"{header} disagrees by {(max_val-max_unc)-(min_val + min_unc)}")
                # print("")


    print("======================Starting combined BAPA figures======================")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    dataset_plot_info = [
        {
            "avg": avg_data_BAPA,
            "err": err_data_BAPA,
            "num": num_data_BAPA,
            "label": "Const. num. particles ($N=300$)",
        },
        {
            "avg": avg_data_CBAPA,
            "err": err_data_CBAPA,
            "num": num_data_CBAPA,
            "label": f"Const. num. fragments ($C={C}$)",
        },
        {
            "avg": avg_data_BAPAWELD,
            "err": err_data_BAPAWELD,
            "num": num_data_BAPAWELD,
            "label": "Const. num. particles, welded ($N=300$)",
        },
    ]

    n_metrics = len(requested_data_headers)
    ncols = 1
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),
        sharex=True,
        constrained_layout=False
    )

    axs = np.asarray(axs).flatten()

    # Same style as your temperature version, but truncate safely
    preferred_ax_order = [2, 3, 0, 1, 4, 5, 6, 7]
    ax_order = [i for i in preferred_ax_order if i < len(axs)]

    # If there are more metrics than entries in preferred_ax_order,
    # append any remaining axes.
    ax_order += [i for i in range(len(axs)) if i not in ax_order]

    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for d_i, info in enumerate(dataset_plot_info):
            avg_data = info["avg"]
            err_data = info["err"]
            num_data = info["num"]

            for t_i, t in enumerate(temps):
                ax.errorbar(
                    M,
                    avg_data[h_i, :, t_i],
                    yerr=err_data[h_i, :, t_i],
                    label=info["label"],
                    color=colors[d_i],       
                    linestyle=styles[d_i],   
                    marker='*' if d_i == 2 else '.',
                    markersize=8,
                    zorder=5
                )

                if include_totals:
                    for m_i, total in enumerate(num_data[h_i, :, t_i]):
                        ax.annotate(
                            f"{total:.0f}",
                            (M[m_i], avg_data[h_i, m_i, t_i]),
                            textcoords="offset points",
                            xytext=(2, 2),
                            fontsize=9,
                            alpha=0.9
                        )

        ax.set_ylabel(label_from_header(header))
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        ax.text(
            0.02, 0.04,
            f"({chr(97 + ax_order[h_i])})",
            transform=ax.transAxes,
            ha='left',
            va='bottom'
        )

    # Common x-label only on bottom row
    for ax in axs[-ncols:]:
        ax.set_xlabel("Fragment size $M$")
    # axs[3].tick_params(labelbottom=True)
    # axs[3].set_xlabel("Fragment size $M$")

    xmin = np.min(M)
    xmax = np.max(M)
    pad = 1.5

    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)

    # Shared legend: get handles from the first used axis
    first_used_ax = axs[ax_order[0]]
    handles, labels = first_used_ax.get_legend_handles_labels()
    labels[0] = labels[0]

    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0001)
    )

    # Remove unused axes
    for i in range(n_metrics, len(axs)):
        fig.delaxes(axs[i])

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    if save_plots:
        plt.savefig(
            f"{figure_folder}{dataset_name}_all_metrics_vs_frag_size.png",
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)  


def add_image_marker(ax, image_path, x, y, zoom=0.12, zorder=10,
                     trim_white=True, white_thresh=245, pad=2):
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    if image_path is None or not os.path.exists(image_path):
        return
    if zoom <= 0:
        return

    if trim_white:
        img = load_and_trim_white(image_path, white_thresh=white_thresh, pad=pad)
    else:
        img = Image.open(image_path).convert("RGB")

    imagebox = OffsetImage(img, zoom=zoom)

    ab = AnnotationBbox(
        imagebox,
        (x, y),
        frameon=False,
        xycoords="data",
        box_alignment=(0.5, 0.5),
        pad=0.0,
        zorder=zorder
    )

    ax.add_artist(ab)


def get_crop_from_white_trim(image_path, white_thresh=245, pad=2):
    """
    Open an image and return crop parameters to take away near-white borders.

    Parameters
    ----------
    image_path : str
        Full path to image.
    white_thresh : int
        Pixels with all RGB channels >= this value are treated as white.
        Lower this if trimming is too aggressive or too weak.
    pad : int
        Number of pixels of padding to add back after cropping.

    Returns
    -------
    img.crop parameters [x0,y0,x1,x2]
    """
    from PIL import Image, ImageOps
    import numpy as np
    import os
    

    img = Image.open(image_path).convert("RGB")
    arr = np.array(img)

    # True where pixel is NOT white
    nonwhite = np.any(arr < white_thresh, axis=2)

    # If image is entirely white, just return original
    if not np.any(nonwhite):
        return img

    ys, xs = np.where(nonwhite)
    x0, x1 = xs.min(), xs.max()
    y0, y1 = ys.min(), ys.max()

    # Add a little padding back
    x0 = max(0, x0 - pad)
    y0 = max(0, y0 - pad)
    x1 = min(arr.shape[1] - 1, x1 + pad)
    y1 = min(arr.shape[0] - 1, y1 + pad)

    return [x0, y0, x1 + 1, y1 + 1]


def load_and_crop(image_path, crop_coords):
    """
    Open an image and crop away given area. Basically wraps img.crop

    Parameters
    ----------
    image_path : str
        Full path to image.
    crop_coords : list

    Returns
    -------
    PIL.Image.Image
        Cropped image.
    """
    from PIL import Image, ImageOps
    img = Image.open(image_path).convert("RGB")
    return img.crop((crop_coords[0], crop_coords[1], crop_coords[2], crop_coords[3]))

def load_and_trim_white(image_path, white_thresh=245, pad=2):
    """
    Open an image and crop away near-white borders.

    Parameters
    ----------
    image_path : str
        Full path to image.
    white_thresh : int
        Pixels with all RGB channels >= this value are treated as white.
        Lower this if trimming is too aggressive or too weak.
    pad : int
        Number of pixels of padding to add back after cropping.

    Returns
    -------
    PIL.Image.Image
        Cropped image.
    """
    from PIL import Image, ImageOps

    crop_coords = get_crop_from_white_trim(image_path,white_thresh,pad)
    
    img = Image.open(image_path).convert("RGB")
    return img.crop((crop_coords[0], crop_coords[1], crop_coords[2], crop_coords[3]))

    
     
def gen_BAPA_plots_images(show_plots=True,save_plots=False,include_totals=False):
    import glob as g
    from PIL import Image
    from matplotlib.lines import Line2D

    image_zoom = 0.02

    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    relax = False

    # job_name = ["DBAPA"]
    job_name = ["BAPA","CBAPA","BAPAWELD","DBAPA"]
    data_prefolders = [path + f'jobs/{name}_' for name in job_name]


    data_prefolder = data_prefolders[0]
    dataset_name = data_prefolder.split("/")[-1]
    figure_folder = path+'data/figures/'

    temps = [1000]
    Nums = [300]
    # M = [1,3,5,10,15,20,30,50,60,75,100,150,300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150,300]
    
    attempts = [i for i in range(30)]


    # requested_data_headers = [1,1,0,1,1,0,0,0,0,1]
    # requested_data_headers = gd.data_headers[:2] + gd.data_headers[3:6] + [gd.data_headers[9]]
    requested_data_headers = [gd.data_headers[0]]
    print(requested_data_headers)
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    
    image_paths = np.full(shape=(len(requested_data_headers),len(M),len(temps)), dtype=object, fill_value="")
    BAPA_image_paths = np.full(shape=image_paths.shape, dtype=object, fill_value="")
    CBAPA_image_paths = np.full(shape=image_paths.shape, dtype=object, fill_value="")
    BAPAWELD_image_paths = np.full(shape=image_paths.shape, dtype=object, fill_value="")
    DBAPA_image_paths = np.full(shape=image_paths.shape, dtype=object, fill_value="")

    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            n = 300
            size = n
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
                            if BAPA_image_paths[h_i,m_i,t_i] == "":
                                blob = f"{path}/data/figures/aggRenders/paper2/ColoredFragg-BAPA_a-*_M-{m}_N-{n}_T-{t}.png"
                                # print(blob)
                                paths = g.glob(blob)
                                # print(paths)
                                if len(paths) == 1:
                                    BAPA_image_paths[h_i,m_i,t_i] = paths[0]

    avg_data_BAPA = np.nanmean(raw_data,axis=1)
    std_data_BAPA = np.nanstd(raw_data,axis=1)
    num_data_BAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPA = std_data_BAPA/np.sqrt(num_data_BAPA)

    # print(avg_data_BAPA)
    
    data_prefolder = data_prefolders[1]
    dataset_name = data_prefolder.split("/")[-1]

    # temps = [1000]
    C = 30
    # # M = [1,3,5,10,15,20,30,50,60,75,100]
    # attempts = [i for i in range(30)]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            if m == 60:
                continue
            n = C*m
            size = n
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            # if h_i == 3:
                            #   print(existing_values_for_size[existing_headers_for_size.index(header)])
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)

                            if CBAPA_image_paths[h_i,m_i,t_i] == "":
                                blob = f"{path}/data/figures/aggRenders/paper2/ColoredFragg-CBAPA_a-*_M-{m}_N-{n}_T-{t}.png"
                                paths = g.glob(blob)
                                if len(paths) == 1:
                                    CBAPA_image_paths[h_i,m_i,t_i] = paths[0]



    avg_data_CBAPA = np.nanmean(raw_data,axis=1)
    std_data_CBAPA = np.nanstd(raw_data,axis=1)
    num_data_CBAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_CBAPA = std_data_CBAPA/np.sqrt(num_data_CBAPA)


    data_prefolder = data_prefolders[2]
    dataset_name = data_prefolder.split("/")[-1]

    # temps = [1000]
    # # M = [3,100]
    # n = 300
    # attempts = [i for i in range(30)]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
                            # raw_data[h_i,a_i,m_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

                            if BAPAWELD_image_paths[h_i,m_i,t_i] == "":
                                blob = f"{path}/data/figures/aggRenders/paper2/ColoredFragg-BAPAWELD_a-*_M-{m}_N-{n}_T-{t}.png"
                                paths = g.glob(blob)
                                if len(paths) == 1:
                                    BAPAWELD_image_paths[h_i,m_i,t_i] = paths[0]

    avg_data_BAPAWELD = np.nanmean(raw_data,axis=1)
    std_data_BAPAWELD = np.nanstd(raw_data,axis=1)
    num_data_BAPAWELD = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPAWELD = std_data_BAPAWELD/np.sqrt(num_data_BAPAWELD)


    C_sp = [1,2,3,4,5,6,10,15,20,30,60] 
    M_sp = 20

    DBAPA_image_paths = np.full(shape=(len(requested_data_headers),len(C_sp),len(temps)), dtype=object, fill_value="")

    data_prefolder = data_prefolders[3]
    dataset_name = data_prefolder.split("/")[-1]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(C_sp),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for c_i,c in enumerate(C_sp):
            n = c*M_sp
            m = M_sp
            size = n
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"

                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            # if h_i == 3:
                            #   print(existing_values_for_size[existing_headers_for_size.index(header)])
                            raw_data[h_i,a_i,c_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)

                            if DBAPA_image_paths[h_i,c_i,t_i] == "":
                                blob = f"{path}/data/figures/aggRenders/paper2/ColoredFragg-DBAPA_a-*_M-{m}_N-{n}_T-{t}.png"
                                paths = g.glob(blob)
                                if len(paths) == 1:
                                    DBAPA_image_paths[h_i,c_i,t_i] = paths[0]



    avg_data_DBAPA = np.nanmean(raw_data,axis=1)
    std_data_DBAPA = np.nanstd(raw_data,axis=1)
    num_data_DBAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_DBAPA = std_data_DBAPA/np.sqrt(num_data_DBAPA)





    print("======================Starting combined BAPA figures======================")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    dataset_plot_info = [
        {
            "avg": avg_data_BAPA,
            "err": err_data_BAPA,
            "num": num_data_BAPA,
            "label": "Const. final size (CFS) $N=300$",
            "image_paths": BAPA_image_paths,
            "xlabel": "Monomers in fragment $M$",
            "xdata": M,
            "legend_label": "$N=300$",
            "legend_loc": [0.5,0.3],
        },
        {
            "avg": avg_data_CBAPA,
            "err": err_data_CBAPA,
            "num": num_data_CBAPA,
            "label": f"Const. num. projectiles (CNP) $C={C}$",
            "image_paths": CBAPA_image_paths,
            "xlabel": "Monomers in fragment $M$",
            "xdata": M,
            "legend_label": "$C=30$",
            "legend_loc": [0.65,0.3],
        },
        {
            "avg": avg_data_BAPAWELD,
            "err": err_data_BAPAWELD,
            "num": num_data_BAPAWELD,
            "label": "Welded const. final size (wCFS) $N=300$",
            "image_paths": BAPAWELD_image_paths,
            "xlabel": "Monomers in fragment $M$",
            "xdata": M,
            "legend_label": "$N=300$",
            "legend_loc": [0.5,0.3],
        },
        {
            "avg": avg_data_DBAPA,
            "err": err_data_DBAPA,
            "num": num_data_DBAPA,
            "label": "DBAPA",
            "image_paths": DBAPA_image_paths,
            "xlabel": "Monomers in final aggregate $N$",
            "xdata": [c*M_sp for c in C_sp],
            "legend_label": "$M=20$",
            "legend_loc": [0.65,0.3],
        },
    ]


    n_metrics = len(requested_data_headers)

    for d_i, info in enumerate(dataset_plot_info):
        avg_data = info["avg"]
        err_data = info["err"]
        num_data = info["num"]
        label = info["label"]
        xlabel = info["xlabel"]
        xdata = info["xdata"]
        legend_label = info["legend_label"]
        legend_loc = info["legend_loc"]

        print(label)
        image_paths = info["image_paths"]


        fig, axs = plt.subplots(
            nrows=1,
            ncols=n_metrics,
            figsize=(7, 5.0),
            sharex=False,
            constrained_layout=False
        )

        axs = np.asarray(axs).flatten()

        ax_order = [i for i in range(len(axs))]

        for h_i, header in enumerate(requested_data_headers):
            ax = axs[ax_order[h_i]]

            for t_i, t in enumerate(temps):

                # Add image at each data point
                for x_i, x in enumerate(xdata):
                    # x = M[m_i]
                    y = avg_data[h_i, x_i, t_i]

                    if image_paths[h_i,x_i,t_i] != "":
                        image_path = image_paths[h_i,x_i,t_i]
                        # print(f"h:{h_i}, m:{m_i}, t:{t_i}")

                        add_image_marker(
                            ax,
                            image_path=image_path,
                            x=x,
                            y=y,
                            zoom=image_zoom,
                            zorder=0
                        )

                        if include_totals:
                            total = num_data[h_i, x_i, t_i]

                            ax.annotate(
                                f"{total:.0f}",
                                (x, y),
                                textcoords="offset points",
                                xytext=(6, 6),
                                fontsize=9,
                                alpha=0.9,
                                zorder=20
                            )
                    else:
                        ax.errorbar(
                            x,
                            avg_data[h_i, x_i, t_i],
                            yerr=err_data[h_i, x_i, t_i],
                            color=colors[t_i],
                            linestyle=None,
                            marker="*",
                            fmt='-',
                            linewidth=1.0,
                            capsize=2,
                            alpha=0.6,
                            zorder=3
                        )

            ax.set_ylabel(label_from_header(header))
            # ax.set_ylabel(header)
            ax.set_xlabel(xlabel)
            ax.set_xscale('log')
            ax.grid(alpha=0.25)

            # ax.text(
            #     0.02, 0.04,
            #     f"({chr(97 + ax_order[h_i])})",
            #     transform=ax.transAxes,
            #     ha='left',
            #     va='bottom'
            # )

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        if len(np.where(~np.isnan(avg_data[0,:,0]))[0]) > 0:
            xmax = xdata[np.where(~np.isnan(avg_data[0,:,0]))[0].max()]
            xmin = xdata[np.where(~np.isnan(avg_data[0,:,0]))[0].min()]
            ymin = np.nanmin(avg_data)
            ymax = np.nanmax(avg_data)
        # xmin = 0
        # xmax = len(M)
        xpad = 1.5
        ypad = 0.03

        for ax in axs[:n_metrics]:
            # ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_xlim(xmin / xpad, xmax * xpad)
            ax.set_ylim(ymin - ypad, ymax + ypad)


        ax.text(
            legend_loc[0], legend_loc[1],
            legend_label,
            fontsize = 30,
            # transform=ax.transAxes,
            transform=fig.transFigure,
            ha='center',
            va='bottom'
        )
        

        if save_plots:
            safe_label = label.replace(" ", "_")
            print(f"{figure_folder}{job_name[d_i]}_image_metrics_vs_frag_size.png")
            plt.savefig(
                f"{figure_folder}{job_name[d_i]}_image_metrics_vs_frag_size.png",
                dpi=300,
                # bbox_inches='tight'
            )

        if show_plots:
            plt.show()

        plt.close(fig)

def permutate_data(data1, data2, permutation=[]):
    assert(len(data1) == len(data2))
    assert(len(data1.shape) == 1)
    assert(len(data2.shape) == 1)
    if len(permutation) == 0:
        permutation = range(len(data1)+len(data2))
    assert(len(data1) + len(data2) == len(permutation))

    all_data = np.concatenate((data1, data2))
    all_data = all_data[permutation]

    data1 = all_data[:len(data1)]
    data2 = all_data[len(data1):len(data1)+len(data2)]

    return data1,data2

    

def test_data(data):
    asym_param_index = 3
    DBAPA_M_ind = 9
    CBAPA_M_ind = 5
    permutations = 50
    # print(asym_param_index)
    # print(gd.data_headers)
    DBAPA_raw_data = data[0,asym_param_index,:,DBAPA_M_ind,0]
    CBAPA_raw_data = data[1,asym_param_index,:,CBAPA_M_ind,0]

    DBAPA_avg_data = np.nanmean(DBAPA_raw_data)
    DBAPA_std_data = np.nanstd(DBAPA_raw_data)
    DBAPA_num_data = np.count_nonzero(~np.isnan(DBAPA_raw_data))
    DBAPA_err_data = DBAPA_std_data/np.sqrt(DBAPA_num_data)
    

    CBAPA_avg_data = np.nanmean(CBAPA_raw_data)
    CBAPA_std_data = np.nanstd(CBAPA_raw_data)
    CBAPA_num_data = np.count_nonzero(~np.isnan(CBAPA_raw_data))
    CBAPA_err_data = CBAPA_std_data/np.sqrt(CBAPA_num_data)
    

    print(f"DBAPA Asymetry Parameter for N=600: {DBAPA_avg_data}+-{DBAPA_err_data}")
    print(f"CBAPA Asymetry Parameter for N=600: {CBAPA_avg_data}+-{CBAPA_err_data}")


    newshape = DBAPA_raw_data.shape

    newshape = (newshape[0]*2)

    all_raw = np.full(shape=(newshape),fill_value=np.nan,dtype=np.float64)
    all_raw[:30] = DBAPA_raw_data
    all_raw[30:] = CBAPA_raw_data

    vals = []
    unc = []
    rng = np.random.default_rng(seed=42)
    for i in range(permutations):
        if i == 0: 
            indices = range(60)
        else:
            indices = rng.permutation(60)

        avg_data_first = np.nanmean(all_raw[indices[0:30]])
        std_data_first = np.nanstd(all_raw[indices[0:30]])
        num_data_first = np.count_nonzero(~np.isnan(all_raw[indices[0:30]]))
        err_data_first = std_data_first/np.sqrt(num_data_first)

        vals.append(avg_data_first)
        unc.append(err_data_first)

        avg_data_sec = np.nanmean(all_raw[indices[30:]])
        std_data_sec = np.nanstd(all_raw[indices[30:]])
        num_data_sec = np.count_nonzero(~np.isnan(all_raw[indices[30:]]))
        err_data_sec = std_data_sec/np.sqrt(num_data_sec)

        # print(f"first Asymetry Parameter for N=600: {avg_data_first}+-{err_data_first}")
        # print(f"secon Asymetry Parameter for N=600: {avg_data_sec}+-{err_data_sec}")

        vals.append(avg_data_sec)
        unc.append(err_data_sec)


    fig, ax = plt.subplots()

    vals = np.asarray(vals)
    unc = np.asarray(unc)

    x = np.arange(len(vals))

    even = x % 2 == 0
    odd = ~even

    ax.errorbar(
        x[even],
        vals[even],
        yerr=unc[even],
        fmt="o",
        linestyle="none",
        color="blue",
    )

    ax.errorbar(
        x[odd],
        vals[odd],
        yerr=unc[odd],
        fmt="o",
        linestyle="none",
        color="red",
    )
    plt.show()

def gen_other_BAPA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    relax = False

    asym_param_index = 3
    DBAPA_M_ind = 9
    CBAPA_M_ind = 5

    data_prefolders = []
    data_prefolders += [path + 'jobs/DBAPA_']
    data_prefolders += [path + 'jobs/CBAPA_']
    figure_folder = path+'data/figures/'


    temps = [1000]
    C_sp = [1,2,3,4,5,6,10,15,20,30,60,61,62]
    M_sp = 20
    C = 30 
    M = [1,3,5,10,15,20,30,50,60,75,100,150,300]

    attempts = [i for i in range(30)]

    requested_data_headers = gd.data_headers[:2] + gd.data_headers[3:6] + [gd.data_headers[9]]
    print(requested_data_headers)

    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    avg_data = np.full(shape=(len(requested_data_headers),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)

    avg_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    std_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)

    avg_data_CBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    std_data_CBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data_CBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data_CBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)

    all_data = np.full(shape=(2,len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)


    for d_i,data_prefolder in enumerate(data_prefolders):
        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
        dataset_name = data_prefolder.split("/")[-1]
        for a_i,a in enumerate(attempts):
            if data_prefolder == data_prefolders[0]: #DBAPA
                xdata = C_sp
            else: #CBAPA
                xdata = M

            for x_i,x in enumerate(xdata):
                if data_prefolder == data_prefolders[0]: #DBAPA
                    m = M_sp
                    c = x
                    n = c*m
                elif data_prefolder == data_prefolders[1]: #CBAPA
                    m = x
                    c = C
                    n = m*c
                size = n
                
                for t_i,t in enumerate(temps):
                    # folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                    # if m == M[0]: #first and last folder is actually from old data
                    #     folder = f"{path}/jobsCosine/lognorm_{a}/N_30/T_{t}/"
                    #     # folder = f"{data_prefolder}{a}/M_{M[0]}/N_{n}/T_{t}/"
                        
                    # else:
                    folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"

                    if os.path.exists(folder+"job_data.csv"):
                        with open(folder+"job_data.csv",'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if n not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(n)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,x_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
                    # elif data_prefolder == data_prefolders[1] and x_i > 6 and x_i < 9:
                    #     print(f"DNE: {folder}job_data.csv")
        # if data_prefolder == == data_prefolders[1]:
        #     np.savetext("rawdata_other.txt",raw_data)
        all_data[d_i] = raw_data
    
    raw_dbapa_permute = all_data[0,asym_param_index,:,DBAPA_M_ind,0]
    raw_cbapa_permute = all_data[1,asym_param_index,:,CBAPA_M_ind,0]
    rng = np.random.default_rng(seed=42)
    permutation = rng.permutation(len(raw_dbapa_permute) + len(raw_cbapa_permute))

    raw_dbapa_permute, raw_cbapa_permute = permutate_data(raw_dbapa_permute, raw_cbapa_permute,permutation)
    
    all_data[0,asym_param_index,:,DBAPA_M_ind,0] = raw_dbapa_permute
    np.savetxt(f"{path}/temp/dbapa_permutation.csv",raw_dbapa_permute)
    all_data[1,asym_param_index,:,CBAPA_M_ind,0] = raw_cbapa_permute
    np.savetxt(f"{path}/temp/cbapa_permutation.csv", raw_cbapa_permute)

    for d_i,data_prefolder in enumerate(data_prefolders):
        avg_data = np.nanmean(all_data[d_i],axis=1)
        std_data = np.nanstd(all_data[d_i],axis=1)
        num_data = np.count_nonzero(~np.isnan(all_data[d_i]),axis=1)
        err_data = std_data/np.sqrt(num_data)

        if data_prefolder == data_prefolders[0]: #DBAPA
            avg_data_DBAPA = avg_data
            std_data_DBAPA = std_data
            num_data_DBAPA = num_data
            err_data_DBAPA = err_data

        elif data_prefolder == data_prefolders[1]: #CBAPA
            avg_data_CBAPA = avg_data
            std_data_CBAPA = std_data
            num_data_CBAPA = num_data
            err_data_CBAPA = err_data

    # #ensure average for repeated sims
    # avg = (avg_data_DBAPA[3,9,0] + avg_data_CBAPA[3,5,0])/2.0
    # unc = np.sqrt(err_data_DBAPA[3,9,0]**2.0 + err_data_CBAPA[3,5,0]**2.0)/2.0

    # avg_data_DBAPA[3,9,0] = avg  
    # err_data_DBAPA[3,9,0] = unc
    # avg_data_CBAPA[3,5,0] = avg
    # err_data_CBAPA[3,5,0] = unc



    # test_data(all_data)
    # exit(0)

    print("======================Starting combined BAPA figures======================")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    dataset_plot_info = [
        {
            "avg": avg_data_DBAPA,
            "err": err_data_DBAPA,
            "num": num_data_DBAPA,
            "label": "Const. frag. size (CFS) $M=20$",
            "color": colors[4],
        },
        {
            "avg": avg_data_CBAPA,
            "err": err_data_CBAPA,
            "num": num_data_CBAPA,
            "label": "Cosnt. num. projectiles (CNP) $C=30$",
            "color": colors[1],
        },
    ]

    n_metrics = len(requested_data_headers)
    ncols = 2
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),
        sharex=True,
        constrained_layout=False
    )

    axs = np.asarray(axs).flatten()

    # Same style as your temperature version, but truncate safely
    preferred_ax_order = [4,3,0,1,2,5]
    preferred_ax_order = [0,1,2,3,4,5]
    preferred_ax_order = [4,5,2,0,1,3]
    preferred_ax_order = [0,1,3,4,5,2]
    ax_order = [i for i in preferred_ax_order if i < len(axs)]
    print([requested_data_headers[i] for i in ax_order])

    # If there are more metrics than entries in preferred_ax_order,
    # append any remaining axes.
    # ax_order += [i for i in range(len(axs)) if i not in ax_order]
        
    xdata_DBAPA = [c*M_sp for c in C_sp]
    xdata_CBAPA = [C*m for m in M]

    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for d_i, info in enumerate(dataset_plot_info):

            if d_i == 0: #DBAPA
                xdata = xdata_DBAPA
            elif d_i == 1: #CBAPA
                xdata = xdata_CBAPA

            avg_data = info["avg"]
            err_data = info["err"]
            num_data = info["num"]
            color = info["color"]

            # print(avg_data[h_i, :, t_i])
            # print(xdata)
            for t_i, t in enumerate(temps):
                ax.errorbar(
                    xdata,
                    avg_data[h_i, :, t_i],
                    yerr=err_data[h_i, :, t_i],
                    label=info["label"],
                    color=color,       
                    linestyle=styles[d_i],   
                    marker='*' if d_i == 2 else '.',
                    markersize=8,
                    zorder=5
                )

                if include_totals:
                    for c_i, total in enumerate(num_data[h_i, :, t_i]):
                        ax.annotate(
                            f"{total:.0f}",
                            (xdata[c_i], avg_data[h_i, c_i, t_i]),
                            textcoords="offset points",
                            xytext=(2, 2),
                            fontsize=9,
                            alpha=0.9
                        )

        ax.set_ylabel(label_from_header(header))
        # ax.set_ylabel(header)
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        ax.text(
            0.02, 0.04,
            f"({chr(97 + ax_order[h_i])})",
            transform=ax.transAxes,
            ha='left',
            va='bottom'
        )

    # Common x-label only on bottom row
    for ax in axs[-ncols:]:
        ax.set_xlabel("Monomers in final aggregate $N$")
    # axs[3].tick_params(labelbottom=True)
    # axs[3].set_xlabel("Fragment size $M$")

    xmin = np.min([np.min(xdata_CBAPA),np.min(xdata_DBAPA)])
    xmax = np.max([np.max(xdata_CBAPA),np.max(xdata_DBAPA)])
    pad = 1.5

    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)

    # Shared legend: get handles from the first used axis
    first_used_ax = axs[ax_order[0]]
    handles, labels = first_used_ax.get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0001)
    )

    # Remove unused axes
    # for i in range(n_metrics, len(axs)):
    #     fig.delaxes(axs[i])

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    if save_plots:
        save_file = f"{figure_folder}all_metrics_vs_N.png"
        print(f"saving to {save_file}")
        plt.savefig(
            save_file,
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)  


def gen_DBAPA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    relax = False

    asym_param_index = 3
    DBAPA_M_ind = 9
    CBAPA_M_ind = 5

    data_prefolders = []
    data_prefolders += [path + 'jobs/DBAPA_']
    data_prefolders += [path + 'jobs/BAPA_']
    figure_folder = path+'data/figures/'


    temps = [1000]
    C_sp = [1,2,3,4,5,6,10,15,20,30,60,61,62]
    M_sp = 20
    C = 30 
    M = [1,3,5,10,15,20,30,50,60,75,100,150,300]
    attempts = [i for i in range(30)]

    requested_data_headers = gd.data_headers[:2] + gd.data_headers[3:6] + [gd.data_headers[9]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    avg_data = np.full(shape=(len(requested_data_headers),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)

    avg_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    std_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data_DBAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)

    avg_data_BAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    std_data_BAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    num_data_BAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)
    err_data_BAPA = np.full(shape=(avg_data.shape),fill_value=np.nan,dtype=np.float64)


    for d_i,data_prefolder in enumerate(data_prefolders):
        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
        dataset_name = data_prefolder.split("/")[-1]
        for a_i,a in enumerate(attempts):
            if data_prefolder == data_prefolders[0]: #DBAPA
                xdata = C_sp
            else: #BAPA
                xdata = M

            for x_i,x in enumerate(xdata):
                if data_prefolder == data_prefolders[0]: #DBAPA
                    m = M_sp
                    c = x
                    n = c*m
                elif data_prefolder == data_prefolders[1]: #BAPA
                    m = x
                    n = 300
                size = n
                
                for t_i,t in enumerate(temps):
                    # folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                    # if m == M[0]: #first and last folder is actually from old data
                    #     folder = f"{path}/jobsCosine/lognorm_{a}/N_30/T_{t}/"
                    #     # folder = f"{data_prefolder}{a}/M_{M[0]}/N_{n}/T_{t}/"
                        
                    # else:
                    folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                    if os.path.exists(folder+"job_data.csv"):
                        with open(folder+"job_data.csv",'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if n not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(n)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                # if h_i == 3:
                                #   print(existing_values_for_size[existing_headers_for_size.index(header)])
                                raw_data[h_i,a_i,x_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)

        if "DBAPA" in data_prefolder: 
            raw_data[asym_param_index,:,DBAPA_M_ind,0] = np.loadtxt(f"{path}/temp/dbapa_permutation.csv")

        avg_data = np.nanmean(raw_data,axis=1)
        std_data = np.nanstd(raw_data,axis=1)
        num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data = std_data/np.sqrt(num_data)

        if data_prefolder == data_prefolders[0]: #DBAPA
            avg_data_DBAPA = avg_data
            std_data_DBAPA = std_data
            num_data_DBAPA = num_data
            err_data_DBAPA = err_data

        elif data_prefolder == data_prefolders[1]: #BAPA
            avg_data_BAPA = avg_data
            std_data_BAPA = std_data
            num_data_BAPA = num_data
            err_data_BAPA = err_data
            


    print("======================Starting combined BAPA figures======================")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    dataset_plot_info = [
        {
            "avg": avg_data_DBAPA,
            "err": err_data_DBAPA,
            "num": num_data_DBAPA,
            "label": "Const. frag. size (CFS) $M=20$",
            "color": colors[4],
        },
        {
            "avg": avg_data_BAPA,
            "err": err_data_BAPA,
            "num": num_data_BAPA,
            "label": "Const. agg. size (CAS) $N=300$",
            "color": colors[0],
        },
    ]

    n_metrics = len(requested_data_headers)
    ncols = 2
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),
        sharex=True,
        constrained_layout=False
    )

    axs = np.asarray(axs).flatten()

    # Same style as your temperature version, but truncate safely
    preferred_ax_order = [4,3,0,1,2,5]
    preferred_ax_order = [0,1,2,3,4,5]
    preferred_ax_order = [4,5,2,0,1,3]
    preferred_ax_order = [0,1,3,4,5,2]
    ax_order = [i for i in preferred_ax_order if i < len(axs)]
    print([requested_data_headers[i] for i in ax_order])

    # If there are more metrics than entries in preferred_ax_order,
    # append any remaining axes.
    # ax_order += [i for i in range(len(axs)) if i not in ax_order]
        
    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for d_i, info in enumerate(dataset_plot_info):

            if d_i == 0: #DBAPA
                xdata = C_sp
            elif d_i == 1: #BAPA
                xdata = [int(300/m) for m in M]

            avg_data = info["avg"]
            err_data = info["err"]
            num_data = info["num"]
            color = info["color"]

            # print(avg_data[h_i, :, t_i])
            # print(xdata)
            for t_i, t in enumerate(temps):
                ax.errorbar(
                    xdata,
                    avg_data[h_i, :, t_i],
                    yerr=err_data[h_i, :, t_i],
                    label=info["label"],
                    color=color,       
                    linestyle=styles[d_i],   
                    marker='*' if d_i == 2 else '.',
                    markersize=8,
                    zorder=5
                )

                if include_totals:
                    for c_i, total in enumerate(num_data[h_i, :, t_i]):
                        ax.annotate(
                            f"{total:.0f}",
                            (xdata[c_i], avg_data[h_i, c_i, t_i]),
                            textcoords="offset points",
                            xytext=(2, 2),
                            fontsize=9,
                            alpha=0.9
                        )

        ax.set_ylabel(label_from_header(header))
        # ax.set_ylabel(header)
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        ax.text(
            0.02, 0.04,
            f"({chr(97 + ax_order[h_i])})",
            transform=ax.transAxes,
            ha='left',
            va='bottom'
        )

    # Common x-label only on bottom row
    for ax in axs[-ncols:]:
        ax.set_xlabel("Number of projectiles $C$")
    # axs[3].tick_params(labelbottom=True)
    # axs[3].set_xlabel("Fragment size $M$")

    xmin = 1#np.min(C)
    xmax = 300#np.max(C)
    pad = 1.5

    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)

    # Shared legend: get handles from the first used axis
    first_used_ax = axs[ax_order[0]]
    handles, labels = first_used_ax.get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0001)
    )

    # Remove unused axes
    # for i in range(n_metrics, len(axs)):
    #     fig.delaxes(axs[i])

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    if save_plots:
        save_file = f"{figure_folder}all_metrics_vs_C.png"
        print(f"saving to {save_file}")
        plt.savefig(
            save_file,
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)  


def gen_BAPA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    relax = False

    asym_param_index = 3
    DBAPA_M_ind = 9
    CBAPA_M_ind = 5

    data_prefolders = [path + 'jobs/BAPA_', path + 'jobs/CBAPA_', path + 'jobs/BAPAWELD_']


    data_prefolder = data_prefolders[0]
    dataset_name = data_prefolder.split("/")[-1]
    figure_folder = path+'data/figures/'

    temps = [1000]
    Nums = [300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150,300]
    
    attempts = [i for i in range(30)]


    # requested_data_headers = [1,1,0,1,1,0,0,0,0,1]
    requested_data_headers = gd.data_headers[:2] + gd.data_headers[3:6] + [gd.data_headers[9]]
    print(requested_data_headers)
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            n = 300
            size = n
            for t_i,t in enumerate(temps):
                # if m == M[-1]:# or m == M[0]: #first and last folder is actually from old data
                #     # folder = f"{path}/jobsCosine/lognorm_{a}/N_300/T_{t}/"
                #     folder = f"{data_prefolder}{a}/M_{M[0]}/N_{n}/T_{t}/"
                    
                # else:
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"

                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)


    avg_data_BAPA = np.nanmean(raw_data,axis=1)
    std_data_BAPA = np.nanstd(raw_data,axis=1)
    num_data_BAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPA = std_data_BAPA/np.sqrt(num_data_BAPA)
    
    # for h_i,header in enumerate(requested_data_headers):
    #     print(header)
    #     print(f"First: {avg_data_BAPA[h_i,0,0]}+-{err_data_BAPA[h_i,0,0]}, Last: {avg_data_BAPA[h_i,-1,0]}+-{err_data_BAPA[h_i,-1,0]}")


    data_prefolder = data_prefolders[1]
    dataset_name = data_prefolder.split("/")[-1]

    temps = [1000]
    C = 30
    # M = [1,3,5,10,15,20,30,50,60,75,100,150]
    attempts = [i for i in range(30)]

    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            n = C*m
            size = n
            # print(f"C: {C}\tm: {m}\tn: {n}")
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            # if h_i == 3:
                            #   print(existing_values_for_size[existing_headers_for_size.index(header)])
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
    if "CBAPA" in data_prefolder: 
        raw_data[asym_param_index,:,CBAPA_M_ind,0] = np.loadtxt(f"{path}/temp/cbapa_permutation.csv")
    # np.savetext("rawdata_other.txt",raw_data)
    avg_data_CBAPA = np.nanmean(raw_data,axis=1)
    std_data_CBAPA = np.nanstd(raw_data,axis=1)
    num_data_CBAPA = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_CBAPA = std_data_CBAPA/np.sqrt(num_data_CBAPA)


    data_prefolder = data_prefolders[2]
    dataset_name = data_prefolder.split("/")[-1]

    temps = [1000]
    # M = [3,100]
    n = 300
    attempts = [i for i in range(30)]

    # requested_data_headers = [1,1,0,1,1,0,0,0,0,1]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[4]]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]

    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,m_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,n,relax)
                            # raw_data[h_i,a_i,m_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

    avg_data_BAPAWELD = np.nanmean(raw_data,axis=1)
    std_data_BAPAWELD = np.nanstd(raw_data,axis=1)
    num_data_BAPAWELD = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data_BAPAWELD = std_data_BAPAWELD/np.sqrt(num_data_BAPAWELD)

    # for m_i,m in enumerate(M):
    #     if m in [3,15,100]:
    #         for h_i,header in enumerate(requested_data_headers):
    #             print(f"BAPA M={m} {header}    : {avg_data_BAPA[h_i,m_i,0]}+-{err_data_BAPA[h_i,m_i,0]} ({num_data_BAPA[h_i,m_i,0]})")
    #             print(f"BAPAWELD M={m} {header}: {avg_data_BAPAWELD[h_i,m_i,0]}+-{err_data_BAPAWELD[h_i,m_i,0]} ({num_data_BAPAWELD[h_i,m_i,0]})")
    #             min_val = np.min([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])
    #             min_unc = [err_data_BAPA[h_i,m_i,0],err_data_BAPAWELD[h_i,m_i,0]][np.argmin([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])]
    #             max_val = np.max([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])
    #             max_unc = [err_data_BAPA[h_i,m_i,0],err_data_BAPAWELD[h_i,m_i,0]][np.argmax([avg_data_BAPA[h_i,m_i,0],avg_data_BAPAWELD[h_i,m_i,0]])]
    #             if min_val + min_unc >= max_val-max_unc:
    #                 print(f"{header} agrees")
    #             else:
    #                 print(f"{header} disagrees by {(max_val-max_unc)-(min_val + min_unc)}")
    #             print("")


    print("======================Starting combined BAPA figures======================")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    dataset_plot_info = [
        {
            "avg": avg_data_BAPA,
            "err": err_data_BAPA,
            "num": num_data_BAPA,
            "label": "Const. agg. size (CAS) $N=300$",
        },
        {
            "avg": avg_data_CBAPA,
            "err": err_data_CBAPA,
            "num": num_data_CBAPA,
            "label": f"Const. num. projectiles (CNP) $C=30$",
        },
        {
            "avg": avg_data_BAPAWELD,
            "err": err_data_BAPAWELD,
            "num": num_data_BAPAWELD,
            "label": "Welded const. agg. size (wCAS) $N=300$",
        },
    ]

    n_metrics = len(requested_data_headers)
    ncols = 2
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),
        sharex=True,
        constrained_layout=False
    )

    axs = np.asarray(axs).flatten()

    # Same style as your temperature version, but truncate safely
    preferred_ax_order = [4,3,0,1,2,5]
    preferred_ax_order = [0,1,2,3,4,5]
    preferred_ax_order = [4,5,2,0,1,3]
    preferred_ax_order = [0,1,3,4,5,2]
    ax_order = [i for i in preferred_ax_order if i < len(axs)]
    print([requested_data_headers[i] for i in ax_order])

    # If there are more metrics than entries in preferred_ax_order,
    # append any remaining axes.
    # ax_order += [i for i in range(len(axs)) if i not in ax_order]

    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for d_i, info in enumerate(dataset_plot_info):
            avg_data = info["avg"]
            err_data = info["err"]
            num_data = info["num"]

            for t_i, t in enumerate(temps):
                ax.errorbar(
                    M,
                    avg_data[h_i, :, t_i],
                    yerr=err_data[h_i, :, t_i],
                    label=info["label"],
                    color=colors[d_i],       
                    linestyle=styles[d_i],   
                    marker='*' if d_i == 2 else '.',
                    markersize=8,
                    zorder=5
                )

                if include_totals:
                    for m_i, total in enumerate(num_data[h_i, :, t_i]):
                        ax.annotate(
                            f"{total:.0f}",
                            (M[m_i], avg_data[h_i, m_i, t_i]),
                            textcoords="offset points",
                            xytext=(2, 2),
                            fontsize=9,
                            alpha=0.9
                        )

        ax.set_ylabel(label_from_header(header))
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        ax.text(
            0.02, 0.04,
            f"({chr(97 + ax_order[h_i])})",
            transform=ax.transAxes,
            ha='left',
            va='bottom'
        )

    # Common x-label only on bottom row
    for ax in axs[-ncols:]:
        ax.set_xlabel("Monomers in projectile $M$")
    # axs[3].tick_params(labelbottom=True)
    # axs[3].set_xlabel("Fragment size $M$")

    xmin = np.min(M)
    xmax = np.max(M)
    pad = 1.5

    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)

    # Shared legend: get handles from the first used axis
    first_used_ax = axs[ax_order[0]]
    handles, labels = first_used_ax.get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0001)
    )

    # Remove unused axes
    # for i in range(n_metrics, len(axs)):
    #     fig.delaxes(axs[i])

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    if save_plots:
        plt.savefig(
            f"{figure_folder}all_metrics_vs_M.png",
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)  
 

def gen_stylized_BAPA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolder = path + 'jobs/BAPA_'

    dataset_name = data_prefolder.split("/")[-1]

    figure_folder = path+'data/figures/'


    temps = [1000]
    # temps = [3,10]
    Nums = [300]
    M = [1,3,5,10,15,20,30,50,60,100]
    
    
    attempts = [i for i in range(30)]

    requested_data_headers = gd.data_headers[:2] + [gd.data_headers[3]] + [gd.data_headers[4]]
    


    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(M),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for m_i,m in enumerate(M):
            for n_i,n in enumerate(Nums):
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
                    if os.path.exists(folder+"job_data.csv"):
                        with open(folder+"job_data.csv",'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if n not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(n)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,m_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)

    avg_data = np.nanmean(raw_data,axis=1)
    std_data = np.nanstd(raw_data,axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data = std_data#/np.sqrt(num_data)

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(avg_data))))
    


    length = len(temps)



    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):
        for n_i,n in enumerate(Nums):
            for t_i,t in enumerate(temps):


                fig, ax = plt.subplots()

                x = M
                mean = avg_data[h_i, :, n_i, t_i]
                sigma = err_data[h_i, :, n_i, t_i]   # if this is 1σ; otherwise replace with your std

                if h_i == 0:
                    color = 'r'
                else:
                    color = colors[h_i]
                plot_with_sigma_bands(
                    ax, x, mean, sigma,
                    label=r'$\mathrm{Avg}$',
                    line_kwargs=dict(linestyle=styles[h_i], marker='.', markersize=10, color=color)
                )

                if include_totals:
                    for txt_i, txt in enumerate(num_data[h_i, :, n_i, t_i]):
                        ax.annotate(f"{txt:0.0f}", (x[txt_i], mean[txt_i]))

                ax.set_xlabel('Fragment size')
                ax.set_ylabel(label_from_header(header))
                # ax.grid(which='major', linewidth=0.6)
                # ax.grid(which='minor', linestyle=':', linewidth=0.5)
                ax.grid(which='major', color='#222222', linewidth=0.6)
                ax.grid(which='minor', color='#222222', linestyle=':', linewidth=0.5)
                    
                ax.set_xscale('log')
                if h_i == 0:
                    ax.legend(loc="lower center")
                elif h_i == 1:
                    ax.legend(loc="lower right")
                elif h_i == 2:
                    ax.legend(loc="lower left")
                # ax.legend()

                plt.tight_layout()

                if save_plots:
                    plt.savefig(f"{figure_folder}{dataset_name}_{header}_stylized_metric_vs_frag_size.png")
                if show_plots:
                    plt.show()
            

def gen_BPCA_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolder = path + 'jobsNovus/constrelax_'
    data_prefolder = path + 'jobsCosine/lognormrelax_'

    dataset_name = data_prefolder.split("/")[-1].strip("_")
    figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

    if save_plots and not os.path.exists(figure_folder):
        os.makedirs(figure_folder)


    temps = [3,10,30,100,300,1000]
    # temps = [3,10]
    Nums = [30,100,300]
    
    
    attempts = [i for i in range(30)]


    # data_shape = (len(M),len(Nums),len(temps))
    


    raw_data = np.full(shape=(len(gd.data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
    for a_i,a in enumerate(attempts):
        for n_i,n in enumerate(Nums):
            for t_i,t in enumerate(temps):
                folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if n not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(n)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(gd.data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)

    avg_data = np.nanmean(raw_data,axis=1)
    std_data = np.nanstd(raw_data,axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data = std_data/np.sqrt(num_data)



    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(avg_data))))
    

    # styles = ['-','--','-.',':']
    # # styles = ['-','--','-.','--.']
    # colors = ['g','b','r','orange','black','red']
    length = len(temps)


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(gd.data_headers):
        for n_i,n in enumerate(Nums):
            for t_i,t in enumerate(temps):

                fig,ax = plt.subplots()


                ax.errorbar(temps,avg_data[h_i,n_i,:],yerr=err_data[h_i,n_i,:],\
                        label=f"N={n},T={t}",color=colors[h_i],\
                        linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

                if include_totals:
                    for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
                        ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

                bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                ax.set_xlabel('Temp (K)')
                ax.set_ylabel(header)
                # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
                # ax.set_xscale('log')
                # if i == 1:
                fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
                plt.tight_layout()
                if save_plots:
                    plt.savefig("{}{}_{}_avgPlot.png".format(figure_folder,dataset_name,header))
                if show_plots:
                    plt.show() 



def gen_BPCA_vs_time_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolder = path + 'jobsNovus/const_'
    data_prefolder = path + 'jobsCosine/lognorm_'

    dataset_name = data_prefolder.split("/")[-1].strip("_")
    figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

    if save_plots and not os.path.exists(figure_folder):
        os.makedirs(figure_folder)


    temp = 3
    # temps = [3,10]
    n = 300

    attempts = [i for i in range(30)]
    
    sizes=list(range(30,301))


    # requested_data_headers = gd.data_headers[:2]
    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]
    requested_data_headers = [gd.data_headers[0]]

    


    for a_i,a in enumerate(attempts):
        raw_data = np.full(shape=(len(requested_data_headers),len(sizes)),fill_value=np.nan,dtype=np.float64)
        folder = f"{data_prefolder}{a}/N_{n}/T_{temp}/"
        print(folder)
        for s_i,size in enumerate(sizes):
            if os.path.exists(folder+"job_data.csv"):
                with open(folder+"job_data.csv",'r') as fp:
                    existing_data = fp.readlines()

                existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                #even though the data can have other sizes in it, 
                #we only want the data of size n
                if size not in existing_sizes:
                    print(f"ERROR: Data of size {n} does not exist for {folder}.")
                    continue
                index = existing_sizes.index(size)*4
                existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                
                for h_i,header in enumerate(requested_data_headers):
                    if header in existing_headers_for_size:
                        raw_data[h_i,s_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)




    
        print("======================Starting figures======================")
        # print(data.shape)
        print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
        



        #   plt.close("all")
        plt.rcParams.update({
            'font.size': 18,
            'text.usetex': True,
            'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
        })

        #Plot metric vs M for all metrics and all N and temps
        for h_i,header in enumerate(requested_data_headers):

            fig,ax = plt.subplots()

            ax.plot(sizes,raw_data[h_i,:],\
                    # label=f"{header} N={n},T={temp}",\
                    color=colors[h_i],\
                    linestyle=styles[h_i],\
                    marker='.',markersize=10,zorder=5)

            # if include_totals:
            #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
            #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            ax.set_xlabel('aggregate size (number of monomers)')



            ax.set_ylabel(label_from_header(header))
            # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
            # ax.set_xscale('log')
            # if i == 1:
            # fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
            plt.tight_layout()
            if save_plots:
                plt.savefig("{}{}_{}_a-{}_t-{}_overtime.png".format(figure_folder,dataset_name,header,a,temp))
            if show_plots:
                plt.show() 
            plt.close()


def gen_BPCA_vs_time_avg_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    # data_prefolder = path + 'jobsCosine/lognorm_'

    # dataset_name = data_prefolder.split("/")[-1]
    # figure_folder = path+'data/figures/BPCA_per_step/'

    data_prefolder = path + 'jobsNovus/const_'
    data_prefolder = path + 'jobsCosine/lognorm_'

    dataset_name = data_prefolder.split("/")[-1].strip("_")
    figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

    if save_plots and not os.path.exists(figure_folder):
        os.makedirs(figure_folder)


    temps = [1000]
    # temps = [3,10]
    n = 300

    attempts = [i for i in range(30)]
    
    sizes=list(range(30,301))


    # requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]
    requested_data_headers = [gd.data_headers[1]]

    


    raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(temps),len(sizes)),fill_value=np.nan,dtype=np.float64)
    print(raw_data.shape)
    for a_i,a in enumerate(attempts):
        for t_i,temp in enumerate(temps):
            folder = f"{data_prefolder}{a}/N_{n}/T_{temp}/"
            # print(folder)
            for s_i,size in enumerate(sizes):
                if os.path.exists(folder+"job_data.csv"):
                    with open(folder+"job_data.csv",'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if size not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(size)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,t_i,s_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                else:
                    pass
                    # print(f"file doesn't exist: {folder}job_data.csv")

    avg_data = np.nanmean(raw_data,axis=1)
    std_data = np.nanstd(raw_data,axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data = std_data/np.sqrt(num_data)

    print(avg_data.shape)
    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
    



    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    for h_i,header in enumerate(requested_data_headers):
        #Plot metric vs M for all metrics and all N and temps
        # for a_i,attempt in enumerate(attempts): 
        fig,ax = plt.subplots()
        for t_i,temp in enumerate(temps):


            ax.errorbar(sizes,avg_data[h_i,t_i,:],yerr=err_data[h_i],\
                    label=f"{header} N={n},T={temp}",color=colors[t_i],\
                    linestyle=styles[h_i],marker='.',markersize=10,zorder=5)


            r0 = 1e-5
            sizes=list(range(30,301))
            Vtot = sizes[0]*(4*np.pi/3)*(r0)**3
            reff = (3*Vtot/(4*np.pi))**(1/3)

            data = Tanaka(sizes,(reff/(1-avg_data[h_i,t_i,0])**(1/3))/(5/3)**(1/2),temp)

            ax.plot(sizes,data,\
                    label=f"Tanaka prediction T={temp}",color=colors[t_i+1],\
                    linestyle=styles[h_i],marker='*',markersize=10,zorder=5)

        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('aggregate size')
        ax.set_ylabel(header)
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        # ax.set_xscale('log')
        # if i == 1:
        fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_avgovertime.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 

    # for h_i,header in enumerate(requested_data_headers):
    #   #Plot metric vs M for all metrics and all N and temps
    #   fig,ax = plt.subplots()
    #   for t_i,temp in enumerate(temps):


    #       ax.errorbar(sizes,avg_data[h_i,t_i,:],yerr=err_data[h_i],\
    #               label=f"{header} N={n},T={temp}",color=colors[t_i],\
    #               linestyle=styles[h_i],marker='.',markersize=10,zorder=5)


    #   # if include_totals:
    #   #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
    #   #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

    #   bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #   ax.set_xlabel('aggregate size')
    #   ax.set_ylabel(header)
    #   # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
    #   # ax.set_xscale('log')
    #   # if i == 1:
    #   fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
    #   plt.tight_layout()
    #   if save_plots:
    #       plt.savefig("{}{}_{}_avgovertime.png".format(figure_folder,dataset_name,header))
    #   if show_plots:
    #       plt.show()


def gen_seqstick_plots(distribution):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    if distribution == "lognormal":
        data_prefolder = path + 'jobs/SeqStickLognormrelax_'
        data_prefolder = path + 'jobsCosine/lognormrelax_'
    elif distribution == "constant":
        data_prefolder = path + 'jobs/SeqStickConstrelax_'
        data_prefolder = path + 'jobsNovus/constrelax_'
    else:
        print("Distribution not recognized")
        exit(-1)

    dataset_name = data_prefolder.split("/")[-1]
    figure_folder = path+'data/figures/'



    n = 300

    attempts = [i for i in range(30)]
    
    size=300

    relax = ("relax" in data_prefolder)
    print(f"relax: {relax}")
    rel = ""
    if relax:
        rel = "relax_"


    requested_data_headers = gd.data_headers[:2] + gd.data_headers[-4:]


    raw_data = np.full(shape=(len(requested_data_headers),len(attempts)),fill_value=np.nan,dtype=np.float64)
    print(f"raw_data shape: {raw_data.shape}")

    for a_i,attempt in enumerate(attempts):
        # folder = f"{data_prefolder}{attempt}/N_{n}/"
        folder = f"{data_prefolder}{attempt}/N_{n}/T_3/"
        # print(folder)
        if os.path.exists(folder+f"{rel}job_data.csv"):
            with open(folder+f"{rel}job_data.csv",'r') as fp:
                existing_data = fp.readlines()

            existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
            #even though the data can have other sizes in it, 
            #we only want the data of size n
            if size not in existing_sizes:
                print(f"ERROR: Data of size {n} does not exist for {folder}.")
                continue
            index = existing_sizes.index(size)*4
            existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
            existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")

            # print(existing_values_for_size)
            
            for h_i,header in enumerate(requested_data_headers):
                if header in existing_headers_for_size:
                    raw_data[h_i,a_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                else:
                    print(f"Header {header} doesnt exist for dir {folder}")
        else:
            print(f"NO DATA FILE FOR FOLDER: {folder}")


    # print(raw_data)
    avg_data = np.nanmean(raw_data,axis=1)
    std_data = np.nanstd(raw_data,axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
    err_data = std_data/np.sqrt(num_data)

    print(f"avg_data shape: {avg_data.shape}")

    for h_i,header in enumerate(requested_data_headers):
        print(f"{distribution} {header}: {avg_data[h_i]} +- {err_data[h_i]} for {num_data[h_i]} data points.")

def gen_BPCA_gcs_csv_tables(save_plots):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolders = []
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')

    for data_prefolder in data_prefolders:
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        table_folder = path+f'data/tables/{dataset_name}/'

        if save_plots and not os.path.exists(table_folder):
            os.makedirs(table_folder)

        temps = [3,10,30,100,300,1000]
        N = [30,100,300]
        attempts = [i for i in range(30)]

        data_file = "job_data.csv" #with centering 

        data_file = "job_data.csv"

        bool_headers = [0,0,0,0,0,0,0,0,0,1]
        requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]
        header = requested_data_headers[0]
        relax = not ("nonrelax" in data_file)
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_path_data_file = folder+f"{rel}{data_file}"
                    if os.path.exists(full_path_data_file):
                        with open(full_path_data_file,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        if header in existing_headers_for_size:
                            raw_data[a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_path_data_file}")

        #convert cm^2 to micrometers^2
        raw_data*=1e8
        avg_data = np.nanmean(raw_data, axis=0)
        std_data = np.nanstd(raw_data, axis=0)

        csv_headers = ["N"]+[f"T{int(T)}K" if float(T).is_integer() else f"T{T}K" for T in temps]
        csv_header_line = ",".join(csv_headers)

        def avg_fmt(x):
            return f"{x:.{3}f}"
        def std_fmt(x):
            return f"{x:.{3}f}"

        filename = os.path.join(table_folder, f"gcsTable-{dataset_name}.csv")
        row = []
        for n_i, n in enumerate(N):
            row.append([f"{n}"]+[f"{avg_fmt(avg_data[n_i, j])} ± {std_fmt(std_data[n_i, j])}" for j in range(len(temps))])

        with open(filename, "w") as f:
            for n_i, n in enumerate(N):
                if n_i == 0:
                    f.write(csv_header_line + "\n")
                f.write(",".join(row[n_i]) + "\n")

        print("Wrote", filename)


def gen_BPCA_porosity_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolders = []
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    data_prefolders.append(path + 'jobsNovus/constrelax_')

    for data_prefolder in data_prefolders:
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        if dataset_name == "constrelax":
            Title = "Constant"
        elif dataset_name == "lognormrelax":
            Title = "Lognormal"
        else:
            Title = ""

        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)


        temps = [3,10,30,100,300,1000]
        N = [300]
        attempts = [i for i in range(30)]


        # data_file = "test_job_data.csv" #without centering #mean mass
        # data_file = "test_maxnc_job_data.csv" #max nc
        # data_file = "DELETE_job_data.csv" #with centering 
        data_file = "job_data.csv" #with centering 
        # data_file = "nonrelax_job_data.csv" 


        data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper
        data_file = "job_data.csv"
        # data_file = "ch32ppb_job_data.csv" 
        # data_file = "ch64ppb_job_data.csv" 
        # data_file = "ch8192ppb_job_data.csv" 




        bool_headers = [0,0,0,0,0,0,0,0,0,1]
        bool_headers = [1,1,0,0,0,0,1,1,1,1]
        # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
        requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

        relax = not ("nonrelax" in data_file)
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_path_data_file = folder+f"{rel}{data_file}"
                    if os.path.exists(full_path_data_file):
                        with open(full_path_data_file,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_path_data_file}")

        avg_data = np.nanmean(raw_data, axis=1)
        std_data = np.nanstd(raw_data, axis=1)
        num_data = np.count_nonzero(~np.isnan(raw_data), axis=1)
        err_data = std_data / np.sqrt(num_data)

        print("======================Starting figures======================")
        # print(data.shape)
        for h_i,header in enumerate(requested_data_headers):
            print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")

        plt.rcParams.update({
            'font.size': 18,
            'text.usetex': True,
            'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
        })

        #Plot metric vs M for all metrics and all N and temps
        fig,ax = plt.subplots(figsize=(10,5))
        color_order = [0,1,4,5,6,7]
        for h_i,header in enumerate(requested_data_headers):

            for n_i,n in enumerate(N):
                # print(avg_data[h_i,n_i,:])
                ax.errorbar(temps,avg_data[h_i,n_i,:],yerr=err_data[h_i,n_i,:],\
                        label=f"{label_from_header(header)}",\
                        color=colors[color_order[h_i]],\
                        linestyle=styles[2],\
                        marker='.',markersize=10,zorder=5)

                if include_totals:
                    for txt_i, txt in enumerate(num_data[h_i,n_i,:]):
                        ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[h_i,n_i,txt_i]))

            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())


        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Porosity')
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        # reserve space for legend
        fig.subplots_adjust(right=0.5)

        handles, labels = ax.get_legend_handles_labels()
        handle_by_label = dict(zip(labels, handles))

        order = np.argsort(avg_data[:,0,-1])[::-1]
        ordered_labels = [label_from_header(requested_data_headers[i]) for i in order]
        ordered_handles = [handle_by_label[l] for l in ordered_labels if l in handle_by_label]


        # put legend in that reserved space (anchored to the axes)
        ax.legend(
            ordered_handles,
            ordered_labels,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),  # 1.02 is just to the right of the axes
            frameon=True,
            borderaxespad=0.0
        )

        h = 0.55
        # start = (0.35, h)    # figure-fraction coords
        end   = (0.68, h)

        fig.text(
            0.45, end[1]+0.03,       # slightly above the head
            Title,
            transform=fig.transFigure,
            color='black',
            fontsize=30,
            ha='center'
        )


        # ---- reserve space for legend ----
        fig.subplots_adjust(right=0.70)

        plt.tight_layout()

        if save_plots:
            plt.savefig(
                f"{figure_folder}{dataset_name}_porosities_overtemp.png",
                dpi=300,
                # bbox_inches='tight'
            )

        if show_plots:
            plt.show()

        plt.close(fig)

def gen_BPCA_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    data_prefolders = []
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    data_prefolders.append(path + 'jobsNovus/constrelax_')

    for data_prefolder in data_prefolders:
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)


        temps = [3,10,30,100,300,1000]
        N = [30,100,300]
        attempts = [i for i in range(30)]

        
        


        # data_file = "test_job_data.csv" #without centering #mean mass
        # data_file = "test_maxnc_job_data.csv" #max nc
        # data_file = "DELETE_job_data.csv" #with centering 
        data_file = "job_data.csv" #with centering 
        # data_file = "nonrelax_job_data.csv" 


        data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper
        data_file = "job_data.csv"
        # data_file = "ch32ppb_job_data.csv" 
        # data_file = "ch64ppb_job_data.csv" 
        # data_file = "ch8192ppb_job_data.csv" 




        bool_headers = [0,0,0,0,0,0,0,0,0,1]
        bool_headers = [1,1,1,1,0,0,1,1,1,1]
        # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
        requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

        relax = not ("nonrelax" in data_file)
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_path_data_file = folder+f"{rel}{data_file}"
                    if os.path.exists(full_path_data_file):
                        with open(full_path_data_file,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)


                    else:
                        print(f"DNE: {full_path_data_file}")


        # print(avg_data)
        plot_ALL_BPCA_vs_temp_plots(raw_data,requested_data_headers,N,temps,figure_folder,dataset_name,include_totals,save_plots,show_plots)
        # plot_individual_BPCA_vs_temp_plots(raw_data,requested_data_headers,N,temps,figure_folder,dataset_name,include_totals)

        # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
        # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")


def plot_ALL_BPCA_vs_temp_plots(
    raw_data,
    requested_data_headers,
    N,
    temps,
    figure_folder,
    dataset_name,
    include_totals=False,
    save_plots=False,
    show_plots=True
):

    avg_data = np.nanmean(raw_data, axis=1)
    std_data = np.nanstd(raw_data, axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data), axis=1)
    err_data = std_data / np.sqrt(num_data)

    print("======================Starting figures======================")
    for h_i, header in enumerate(requested_data_headers):
        print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    # ---- layout ----
    n_metrics = len(requested_data_headers)
    ncols = 2
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),      # page-sized
        sharex=True,
        constrained_layout=True
    )

    axs = axs.flatten()

    # ---- plotting ----
    ax_order = [2,3,0,1,4,5,6,7]
    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for n_i, n in enumerate(N):
            ax.errorbar(
                temps,
                avg_data[h_i, n_i, :],
                yerr=err_data[h_i, n_i, :],
                label=f"N={n}",
                color=colors[h_i],      # assumes defined elsewhere
                linestyle=styles[n_i],  # assumes defined elsewhere
                marker='.',
                markersize=8,
                zorder=5
            )

            if include_totals:
                for t_i, total in enumerate(num_data[h_i, n_i, :]):
                    ax.annotate(
                        f"{total:.0f}",
                        (temps[t_i], avg_data[h_i, n_i, t_i]),
                        textcoords="offset points",
                        xytext=(2, 2),
                        fontsize=9,
                        alpha=0.9
                    )

        ax.set_ylabel(label_from_header(header))
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        # panel label (a), (b), ...
        ax.text(
            0.02, 0.04,
            f"({chr(97 + ax_order[h_i])})",
            transform=ax.transAxes,
            ha='left', va='bottom'
        )

    # ---- common x-label only on bottom row ----
    for ax in axs[-ncols:]:
        ax.set_xlabel('Temperature [K]')

    xmin = np.min(temps)
    xmax = np.max(temps)

    pad = 1.5
    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)



    # ---- single shared legend ----
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=len(N),
        frameon=False,
        bbox_to_anchor=(0.5,1.0001)
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])


    # ---- remove unused axes ----
    for i in range(n_metrics, len(axs)):
        fig.delaxes(axs[i])

    # ---- save / show ----
    if save_plots:
        plt.savefig(
            f"{figure_folder}{dataset_name}_all_metrics_overtemp.png",
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)

def plot_individual_BPCA_vs_temp_plots(
    raw_data,
    requested_data_headers,
    N,
    temps,
    figure_folder,
    dataset_name,
    include_totals=False,
    save_plots=False,
    show_plots=True
):

    avg_data = np.nanmean(raw_data, axis=1)
    std_data = np.nanstd(raw_data, axis=1)
    num_data = np.count_nonzero(~np.isnan(raw_data), axis=1)
    err_data = std_data / np.sqrt(num_data)

    print("======================Starting figures======================")
    # print(data.shape)
    for h_i,header in enumerate(requested_data_headers):
        print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")
    

    # styles = ['-','--','-.',':']
    # # styles = ['-','--','-.','--.']
    # colors = ['g','b','r','orange','black','red']


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])
            ax.errorbar(temps,avg_data[h_i,n_i,:],yerr=err_data[h_i,n_i,:],\
                    label=f"N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)

            if include_totals:
                for txt_i, txt in enumerate(num_data[h_i,n_i,:]):
                    ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[h_i,n_i,txt_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')


        # print(header)
        # print(label_from_header(header))
        ax.set_ylabel(label_from_header(header))
        # ax.set_title(f'Constant size distribution asymmetry vs temp')
        ax.set_xscale('log')
        # if header == requested_data_headers[1]:
        #   fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_overtemp.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()


def gen_BPCA_ratio_bugbetter_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    N = [30,100,300]
    attempts = [i for i in range(30)]

    
    data_files = []
    data_files.append("job_data.csv")
    data_files.append("test_job_data.csv") #without centering #mean mass



    bool_headers = [1,1,0,0]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    # data_prefolders = []
    # data_prefolders.append(path + 'jobsNovus/constrelax_')
    # data_prefolders.append(path + 'jobsCosine/lognormrelax_')

    data_prefolder = path + 'jobsNovus/constrelax_'
    # data_prefolder = path + 'jobsCosine/lognormrelax_'
    
    avg_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    for d_i,data_file in enumerate(data_files):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)

        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_data_path = folder+f"{rel}{data_file}"
                    if os.path.exists(full_data_path):
                        print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_data_path}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


    ratio_data = avg_data[0]/avg_data[1]
    ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


    # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
    # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
    

    


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])
            # ax.plot(temps,ratio_data[h_i,n_i,:],\
            ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
                    label=f"N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)

        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')

        ax.axhline(1)

        ax.set_ylabel(label_from_header(header))
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        ax.set_xscale('log')
        if header == requested_data_headers[-1]:
            fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_ratiobettertobugovertemp.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()



def gen_BPCA_ratio_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    N = [30,100,300]
    attempts = [i for i in range(30)]

    

    data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"
    # data_file = "ch64ppb_job_data.csv" 



    bool_headers = [1,1,1,1,0,0,1,1,1,1]
    bool_headers = [1,1,1,1,0,0,1,1,1,1]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolders = []
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    
    avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = not ("nonrelax" in data_file)

        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_data_path = folder+f"{rel}{data_file}"
                    if os.path.exists(full_data_path):
                        print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size: 
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                                

                                # pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                                # r = np.power((1/size)*np.sum(radius**3),1.0/3.0)
                                # raw_data[h_i,a_i,n_i,t_i] = float(existing_values_for_size[existing_headers_for_size.index(header)])*size**(1.0/3.0)
                                
                                # pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                                # # r = np.power((1/size)*np.sum(radius**3),1.0/3.0)
                                # S = float(existing_values_for_size[existing_headers_for_size.index(header)])*((np.pi*np.sum(np.power(radius,2))))
                                # r = np.sqrt(S/np.pi)
                                # # print(f"r: {r}")
                                # r_ef_cubed = np.sum(np.power(radius,3))
                                # # print(f"r_ef: {r_ef_cubed**(1.0/3.0)}")
                                # data = 1-(r_ef_cubed/r**3)
                                # raw_data[h_i,a_i,n_i,t_i] = data
                                
                                # pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                                # S = float(existing_values_for_size[existing_headers_for_size.index(header)])*((np.pi*np.sum(np.power(radius,2))))
                                # r = np.power((1/size)*np.sum(radius**3),1.0/3.0)
                                # data = S*size**(1.0/3.0)/(size*np.pi*r**2)
                                # raw_data[h_i,a_i,n_i,t_i] = data

                    else:
                        print(f"DNE: {full_data_path}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


    ratio_data = avg_data[0]/avg_data[1]
    ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)

    for h_i, header in enumerate(requested_data_headers):
        print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")

    plot_ALL_BPCA_ratio_vs_temp_plots(ratio_data,ratio_errs,requested_data_headers,N,temps,figure_folder,dataset_name,include_totals,save_plots,show_plots)
    # plot_individual_BPCA_ratio_vs_temp_plots(ratio_data,ratio_errs,requested_data_headers,N,temps,figure_folder,dataset_name,include_totals,save_plots,show_plots)


    # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
    # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

def plot_ALL_BPCA_ratio_vs_temp_plots(
    ratio_data,
    ratio_errs,
    requested_data_headers,
    N,
    temps,
    figure_folder,
    dataset_name,
    include_totals=False,
    save_plots=False,
    show_plots=True
):


    print("======================Starting figures======================")
    

    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    # ---- layout ----
    n_metrics = len(requested_data_headers)
    ncols = 2
    nrows = math.ceil(n_metrics / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(8.5, 11.0),      # page-sized
        sharex=True,
        constrained_layout=True
    )

    axs = axs.flatten()

    labels = [label_from_header(i) for i in requested_data_headers]

    # ---- plotting ----
    ax_order = [2,3,0,1,4,5,6,7]
    for h_i, header in enumerate(requested_data_headers):
        ax = axs[ax_order[h_i]]

        for n_i, n in enumerate(N):
            ax.errorbar(
                temps,
                ratio_data[h_i, n_i, :],
                yerr=ratio_errs[h_i, n_i, :],
                label=f"N={n}",
                color=colors[h_i],      # assumes defined elsewhere
                linestyle=styles[n_i],  # assumes defined elsewhere
                marker='.',
                markersize=8,
                zorder=5
            )
            ax.axhline(1)

            ax.text(0.97, 0.95, labels[h_i],
                transform=ax.transAxes,  # now (0,0) = bottom-left, (1,1) = top-right of the axes
                ha="right", va="top")

            if include_totals:
                for t_i, total in enumerate(num_data[h_i, n_i, :]):
                    ax.annotate(
                        f"{total:.0f}",
                        (temps[t_i], avg_data[h_i, n_i, t_i]),
                        textcoords="offset points",
                        xytext=(2, 2),
                        fontsize=9,
                        alpha=0.9
                    )

        ax.set_ylabel("Ratio")
        ax.set_xscale('log')
        ax.grid(alpha=0.25)

        # panel label (a), (b), ...
        # ax.text(
        #   0.02, 0.04,
        #   f"({chr(97 + ax_order[h_i])})",
        #   transform=ax.transAxes,
        #   ha='left', va='bottom'
        # )

    # ---- common x-label only on bottom row ----
    for ax in axs[-ncols:]:
        ax.set_xlabel('Temperature [K]')

    xmin = np.min(temps)
    xmax = np.max(temps)

    pad = 1.5
    for ax in axs[:n_metrics]:
        ax.set_xlim(xmin / pad, xmax * pad)



    # ---- single shared legend ----
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc='upper center',
        ncol=len(N),
        frameon=False,
        bbox_to_anchor=(0.5,1.0001)
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])


    # ---- remove unused axes ----
    for i in range(n_metrics, len(axs)):
        fig.delaxes(axs[i])

    # ---- save / show ----
    if save_plots:
        plt.savefig(
            f"{figure_folder}all_metrics_ratio_overtemp.png",
            dpi=300,
            bbox_inches='tight'
        )

    if show_plots:
        plt.show()

    plt.close(fig)


def plot_individual_BPCA_ratio_vs_temp_plots(
    ratio_data,
    ratio_errs,
    requested_data_headers,
    N,
    temps,
    figure_folder,
    dataset_name,
    include_totals,
    save_plots,show_plots):

    print("======================Starting figures======================")

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    labels = [label_from_header(i) for i in requested_data_headers]
    
    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])
            # ax.plot(temps,ratio_data[h_i,n_i,:],\
            ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
                    label=f"N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)

        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')

        ax.axhline(1)

        ax.text(0.82, 0.95, labels[h_i],
            transform=ax.transAxes,  # now (0,0) = bottom-left, (1,1) = top-right of the axes
            ha="left", va="top")
        # if h_i == 0:
        #   ax.text(400,1.19,labels[h_i])
        # elif h_i == 1:
        #   ax.text(400,1.1,labels[h_i])
        # elif h_i == 2:
        #   ax.text(400,1.035,labels[h_i])
        # elif h_i == 3:
        #   ax.text(400,1.01,labels[h_i])


        # ax.set_ylabel(label_from_header(header))
        ax.set_ylabel("Ratio")
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        ax.set_xscale('log')
        # if header == requested_data_headers[-1]:
            # fig.legend(loc='lower right',bbox_to_anchor=(0.97, 0.96))
        fig.legend(loc='lower right',bbox_to_anchor=(0.97, 0.165))

        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_ratioovertemp.png".format(figure_folder,header))
    if show_plots:
        plt.show() 
    plt.close()

def gen_BPCA_double_ratio_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    N = [30,100,300]
    attempts = [i for i in range(30)]

    
    data_files = []
    data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_files.append("job_data.csv")



    bool_headers = [1,1,0,1]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolders = []
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    
    ratio_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    ratio_errs = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    double_ratio_data = np.full(shape=(len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    for df_i,data_file in enumerate(data_files):
        avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        std_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        num_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        err_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for d_i,data_prefolder in enumerate(data_prefolders):
            dataset_name = data_prefolder.split("/")[-1].strip("_")
            figure_folder = path+f'data/figures/{dataset_name}/'

            if save_plots and not os.path.exists(figure_folder):
                os.makedirs(figure_folder)

            relax = False
            if data_file == "job_data.csv":
                relax = True
            print(f"relax: {relax}")
            rel = ""
            if relax:
                rel = "relax_"

            raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
            for a_i,a in enumerate(attempts):
                for n_i,n in enumerate(N):
                    size = n
                    for t_i,t in enumerate(temps):
                        folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                        full_data_path = folder+f"{rel}{data_file}"
                        if os.path.exists(full_data_path):
                            print(f"opening {full_data_path}")
                            with open(full_data_path,'r') as fp:
                                existing_data = fp.readlines()

                            existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                            #even though the data can have other sizes in it, 
                            #we only want the data of size n
                            if size not in existing_sizes:
                                print(f"ERROR: Data of size {n} does not exist for {folder}.")
                                continue
                            index = existing_sizes.index(size)*4
                            existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                            existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                            
                            for h_i,header in enumerate(requested_data_headers):
                                if header in existing_headers_for_size:
                                    raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                        else:
                            print(f"DNE: {full_data_path}")


            avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
            std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
            num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
            err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


        ratio_data[df_i] = avg_data[0]/avg_data[1]
        ratio_errs[df_i] = ratio_data[df_i]*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


    double_ratio_data = ratio_data[0]/ratio_data[1]
    double_ratio_errs = double_ratio_data*np.sqrt((ratio_errs[0]/ratio_data[0])**2+(ratio_errs[1]/ratio_data[1])**2)
    # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
    # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
    


    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])

            # ax.plot(temps,double_ratio_data[h_i,n_i,:],\
            ax.errorbar(temps,double_ratio_data[h_i,n_i,:],yerr=double_ratio_errs[h_i,n_i,:],\
                    label=f"N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)

        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')

        ax.axhline(1)

        ax.set_ylabel(label_from_header(header))
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        ax.set_xscale('log')
        if header == requested_data_headers[-1]:
            fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_doubleratioovertemp.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()



def gen_BPCA_ratio_nonreltorel_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    N = [30,100,300]
    attempts = [i for i in range(30)]

    
    data_files = []
    data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_files.append("job_data.csv")


    bool_headers = [1,1,1,1]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolder = path + 'jobsNovus/constrelax_'
    data_prefolder = path + 'jobsCosine/lognormrelax_'
    
    avg_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    for d_i,data_file in enumerate(data_files):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_data_path = folder+f"{rel}{data_file}"
                    if os.path.exists(full_data_path):
                        print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_data_path}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


    ratio_data = avg_data[0]/avg_data[1]
    ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


    # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
    # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])
            # ax.plot(temps,ratio_data[h_i,n_i,:],\
            ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
                    label=f"N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)

        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')

        ax.axhline(1)

        ax.set_ylabel(label_from_header(header))
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        ax.set_xscale('log')
        if header == requested_data_headers[-1]:
            fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_rationonreltorelovertemp.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()

def gen_BPCA_temp_sensitivity_plots(show_plots=True,save_plots=False,include_totals=False):
    def S(sigma):
        return np.sum(1/sigma**2)

    def Si(i,sigma):
        return np.sum(i/(sigma**2))

    def Sii(x,y,sigma):
        return np.sum(x*y/(sigma**2))

    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    x = np.log(temps)
    N = [30,100,300]
    attempts = [i for i in range(30)]

    
    data_file = "test_RKBMs_job_data.csv" #print both ways of RKBM
    data_file = "test_job_data.csv" #without centering #mean mass
    data_file = "test_maxnc_job_data.csv" #max nc
    data_file = "job_data.csv" #with centering
    # data_file = "ch64ppb_job_data.csv" 



    bool_headers = [0,0,0,0,0,0,0,0,0,1]
    bool_headers = [1,1,1,1,0,0,1,1,1,1]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    order = [2,3,0,1,4,5,6,7]
    requested_data_headers = [requested_data_headers[i] for i in order]

    data_prefolders = []
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    
    for data_prefolder in data_prefolders:
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        if dataset_name == "constrelax":
            Title = "Constant"
        elif dataset_name == "lognormrelax":
            Title = "Lognormal"
        else:
            Title = ""

        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = not ("nonrelax" in data_file)
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        slope_data = np.full(shape=(len(requested_data_headers),len(N)),fill_value=0,dtype=np.float64)
        slope_sigma_data = np.full(shape=(len(requested_data_headers),len(N)),fill_value=0,dtype=np.float64)
        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    if os.path.exists(folder+f"{rel}{data_file}"):
                        full_data_path = folder+f"{rel}{data_file}"
                        # print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                                

                                # raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

                                # pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                                # data = float(existing_values_for_size[existing_headers_for_size.index(header)])*(np.pi*np.sum(np.power(radius,2)))/size**(1.0/3.0)
                                # raw_data[h_i,a_i,n_i,t_i] = data
                                
                                # raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]


        avg_data = np.nanmean(raw_data,axis=1)
        std_data = np.nanstd(raw_data,axis=1)
        num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data = std_data/np.sqrt(num_data)

        for h_i,header in enumerate(requested_data_headers):
            for n_i,n in enumerate(N):
                y = avg_data[h_i,n_i,:]
                if (h_i == requested_data_headers[-1]):
                    y = [n**(1.0/3.0)*i for i in y]
                sigma = err_data[h_i,n_i,:]

                delta = S(sigma)*Sii(x,x,sigma)-(Si(x,sigma))**2
                slope_data[h_i,n_i] = np.abs((S(sigma)*Sii(x,y,sigma)-Si(x,sigma)*Si(y,sigma))/delta)
                slope_sigma_data[h_i,n_i] = np.sqrt(S(sigma)/delta)
        # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
        # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

        
        print("======================Starting figures======================")
        # print(data.shape)
        print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
        


        nMetrics = len(requested_data_headers)
        nN = len(N)

        plt.rcParams.update({
            'font.size': 18,
            'text.usetex': True,
            'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
        })

        # fig, ax = plt.subplots()
        fig, ax = plt.subplots(figsize=(15,5))

        # ----------------------------------------------------
        # Automatic x-position generation
        # ----------------------------------------------------
        panel_width = 1.0 / nMetrics
        panel_centers = (np.arange(nMetrics) + 0.5) * panel_width
        jitter = np.linspace(-panel_width*2/9, panel_width*2/9, nN)   # horizontal offsets for each N value

        # Store tick info
        xticks = []
        xticklabels = []

        # ----------------------------------------------------
        # Plot each metric (panel)
        # ----------------------------------------------------
        for m in range(nMetrics):

            center = panel_centers[m]

            for i, n in enumerate(N):

                # Compute x-locations for this metric's points (with jitter)
                x_plot = center + jitter[i]

                ax.errorbar(
                    x_plot,
                    slope_data[m, i],
                    yerr=slope_sigma_data[m, i],
                    fmt='o',
                    linewidth=2,
                    capsize=6,
                    color=colors[i]
                )

                # Save ticks and labels
                xticks.append(x_plot)
                xticklabels.append(str(n))

            # Add vertical divider between panels (except the first)
            if m > 0:
                ax.axvline(
                    x=center - panel_width/2,
                    color='black'
                )

            ax.text(
                panel_centers[m],                   # x in axes coordinates
                0.95,                     # y near the top
                label_from_header(requested_data_headers[m]),
                # fontsize=22,
                ha='center', va='top',
                transform=ax.transAxes    # <-- IMPORTANT
            )

        # ----------------------------------------------------
        # Zero horizontal line
        # ----------------------------------------------------
        ax.axhline(0, linestyle='--', color='black')

        # ----------------------------------------------------
        # Set ticks, labels, limits
        # ----------------------------------------------------
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)

        ax.set_xlim(0,1)
        # ax.set_xlim(0 - panel_width, 1 + panel_width)

        # ax.set_ylim(
        #   np.min(slope_data - slope_sigma_data) - 0.002,
        #   np.max(slope_data + slope_sigma_data) + 0.005
        # )
        ax.set_ylim(-0.0010770496039509136, 0.02782241381226227)

        ax.set_xlabel("Aggregate size (N)")
        ax.set_ylabel("Sensitivity to Temperature")

        fig.tight_layout()


        h = 0.61
        start = (0.35, h)    # figure-fraction coords
        end   = (0.68, h)

        # arrow = FancyArrowPatch(
        #     start, end,
        #     transform=fig.transFigure,
        #     arrowstyle='-|>',
        #     mutation_scale=25,
        #     lw=5, color='black',
        #     zorder=200
        # )
        # fig.add_artist(arrow)

        fig.text(
            # 0.515, end[1]+0.03,       # slightly above the head
            0.7, end[1]+0.03,       # slightly above the head
            Title,
            transform=fig.transFigure,
            color='black',
            fontsize=30,
            ha='center'
        )

        if save_plots:
            plt.savefig("{}{}_methodComp.png".format(figure_folder,dataset_name))
        if show_plots:
            plt.show() 
        plt.close()


def gen_BPCA_rolling_fric_plots(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    N = [300]
    attempts = [i for i in range(30)]

    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    bool_headers = [1,1,0,0,1,1]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolders = []
    data_prefolders.append(path + 'jobs/constrollingfricrelax')
    # data_prefolders.append(path + 'jobs/constrollingfric')
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    
    avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_data_path = folder+f"{rel}{data_file}"
                    if os.path.exists(full_data_path):
                        print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_data_path}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


    # ratio_data = avg_data[0]/avg_data[1]
    # ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


    # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
    # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

    
    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    for h_i,header in enumerate(requested_data_headers):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[h_i,n_i,:])
            # ax.plot(temps,ratio_data[h_i,n_i,:],\
            ax.errorbar(temps,avg_data[0,h_i,n_i,:],yerr=err_data[0,h_i,n_i,:],\
                    label=f"rolling fric 0.001 N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[0],\
                    marker='.',markersize=10,zorder=5)
            if include_totals:
                for txt_i, txt in enumerate(num_data[0,h_i,n_i,:]):
                    ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[0,h_i,n_i,txt_i]))

            ax.errorbar(temps,avg_data[1,h_i,n_i,:],yerr=err_data[1,h_i,n_i,:],\
                    # label=f"rolling fric 0.001 not relaxed N={n}",\
                    label=f"rolling fric 1e-5 N={n}",\
                    color=colors[h_i],\
                    linestyle=styles[1],\
                    marker='.',markersize=10,zorder=5)

            if include_totals:
                for txt_i, txt in enumerate(num_data[1,h_i,n_i,:]):
                    ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[1,h_i,n_i,txt_i]))


        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_xlabel('Temperature in K')

        ax.set_ylabel(label_from_header(header))
        # ax.set_title('{} {} vs Temp'.format(dataset_name,method))
        # ax.set_title('Both relaxed')
        ax.set_xscale('log')
        # if header == requested_data_headers[-1]:
        fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_rationonreltorelovertemp.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()



def gen_BPCA_porosity_vs_asymmetry(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    temps = [3,10,30,100,300,1000]
    attempts = [i for i in range(30)]
    N = [30,100,300]

    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    bool_headers = [1,0,0,0,1,0]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolders = []
    # data_prefolders.append(path + 'jobs/constrollingfricrelax')
    # data_prefolders.append(path + 'jobs/constrollingfric')
    data_prefolders.append(path + 'jobsNovus/constrelax_')
    data_prefolders.append(path + 'jobsCosine/lognormrelax_')
    # data_prefolders.append(path + 'jobsNovus/const_')
    # data_prefolders.append(path + 'jobsCosine/lognorm_')
    
    avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            for n_i,n in enumerate(N):
                size = n
                for t_i,t in enumerate(temps):
                    folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
                    full_data_path = folder+f"{rel}{data_file}"
                    if os.path.exists(full_data_path):
                        print(f"opening {full_data_path}")
                        with open(full_data_path,'r') as fp:
                            existing_data = fp.readlines()

                        existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                        #even though the data can have other sizes in it, 
                        #we only want the data of size n
                        if size not in existing_sizes:
                            print(f"ERROR: Data of size {n} does not exist for {folder}.")
                            continue
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                        
                        for h_i,header in enumerate(requested_data_headers):
                            if header in existing_headers_for_size:
                                raw_data[h_i,a_i,n_i,t_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                    else:
                        print(f"DNE: {full_data_path}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


        # ratio_data = avg_data[0]/avg_data[1]
        # ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


        # print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
        # print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")


        print("======================Starting figures======================")
        # print(data.shape)
        print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

        #   plt.close("all")
        plt.rcParams.update({
            'font.size': 18,
            'text.usetex': True,
            'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
        })

        #Plot metric vs M for all metrics and all N and temps
        # for t_i,temp in enumerate(temps):

        fig,ax = plt.subplots()

        for n_i,n in enumerate(N):
            # print(avg_data[d_i,0,n_i,:])
            # print(avg_data[d_i,1,n_i,:])
            # ax.plot(temps,ratio_data[h_i,n_i,:],\
            ax.errorbar(avg_data[d_i,1,n_i,:],avg_data[d_i,0,n_i,:],\
                    yerr=err_data[d_i,0,n_i,:],\
                    xerr=err_data[d_i,1,n_i,:],\
                    label=f"N={n}",\
                    color=colors[d_i],\
                    linestyle=styles[n_i],\
                    marker='.',markersize=10,zorder=5)
            if include_totals:
                for txt_i, txt in enumerate(num_data[0,h_i,n_i,:]):
                    ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[0,h_i,n_i,txt_i]))



        # if include_totals:
        #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
        #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        ax.set_ylabel(label_from_header(requested_data_headers[0]))
        ax.set_xlabel(label_from_header(requested_data_headers[1]))

        # ax.set_ylabel(label_from_header(header))
        # ax.set_title(f'{data_prefolder}')
        # ax.set_title('Both relaxed')
        # ax.set_xscale('log')
        # if header == requested_data_headers[-1]:
        fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))

        #Shade the region:
        # Save current limits so autoscaling from fill doesn't move them
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        # Domain must be A>0 (avoid division by zero)
        if ax.get_xscale() == 'log':
            A = np.geomspace(max(xmin, 1e-12), xmax, 500)
        else:
            A = np.linspace(max(xmin, 1e-12), xmax, 500)


        Aminus_one = A-np.full_like(A,1)

        P_boundry = 1-np.power(0.49/Aminus_one,3.0/4.0)
        # print(xmin,xmax)
        # print(P_boundry[0],P_boundry[-1])
        # print(1-np.power(0.49/(12-1),3.0/4.0))
        # print(1+0.49/np.power(1-0.4,4.0/3.0))

        # Draw boundary
        ax.plot(A, P_boundry, linestyle='--', linewidth=1.5, zorder=3)
        # ax.plot(A, P_boundry, linestyle='--', linewidth=1.5, zorder=3, label=r"$\mathcal{P}_{abc}=1-\left(\frac{0.49}{\mathcal{A}-1}\right)^{3/4}$")

        # Shade region A > 1 + 2/P (everything above the curve)
        ax.fill_between(A, P_boundry, 0, alpha=0.12, zorder=0, label="Allowable region")

        # Restore limits so fill doesn't change view
        # ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, ymax)


        plt.tight_layout()
        if save_plots:
            plt.savefig("{}{}_{}_porosity_vs_asymmetry.png".format(figure_folder,dataset_name,header))
        if show_plots:
            plt.show() 
        plt.close()


def gen_c_over_a_plot(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    attempts = [i for i in range(30)]
    N = [300]
    M = [1,3,5,10,15,20,30,50,60,100,150]
    C=30
    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    bool_headers = [1,0,0,0,1,0,0,0,0,0]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = ["c over a"]

    data_prefolders = []

    data_prefolders.append(path + 'jobs/BAPA_')
    data_prefolders.append(path + 'jobs/BAPAWELD_')
    data_prefolders.append(path + 'jobs/CBAPA_')

    avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    BAPA_avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    BAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    BAPAWELD_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    CBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
        for a_i,att in enumerate(attempts):
            size = N[0]
            n_i = 0
            n=size

            
            for m_i,m in enumerate(M):
                h_i = 0
                if data_prefolder == data_prefolders[2]: #CBAPA
                    size = C*m
                    n = size
                folder = f"{data_prefolder}{att}/M_{m}/N_{n}/T_1000/"
                print(folder)
                # full_data_path = folder+f"{rel}{data_file}"
                if os.path.exists(f'{folder}{n}_simData.csv'):
                    pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                    if pos is None:
                        raw_data[h_i,a_i,n_i,m_i] = np.nan
                    else:
                        effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 
                        a,b,c = u.calc_equivalent_ellipsoid_principal_axes(effective_radius,pos,mass,folder,write=False)
                        raw_data[h_i,a_i,n_i,m_i] = c/a
        
                else:
                    print(f"DNE: {folder}")


        avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])

        if data_prefolder == data_prefolders[0]:
            BAPA_avg_data = avg_data
            BAPA_std_data = std_data
            BAPA_num_data = num_data
            BAPA_err_data = err_data
        elif data_prefolder == data_prefolders[1]:
            BAPAWELD_avg_data = avg_data
            BAPAWELD_std_data = std_data
            BAPAWELD_num_data = num_data
            BAPAWELD_err_data = err_data
        elif data_prefolder == data_prefolders[2]:
            CBAPA_avg_data = avg_data
            CBAPA_std_data = std_data
            CBAPA_num_data = num_data
            CBAPA_err_data = err_data


    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    # for t_i,temp in enumerate(temps):

    fig,ax = plt.subplots()

    for d_i,data_prefolder in enumerate(data_prefolders):
        if data_prefolder == data_prefolders[0]:
            avg_data = BAPA_avg_data
            std_data = BAPA_std_data 
            num_data = BAPA_num_data 
            err_data = BAPA_err_data 
            label = "BAPA"
            color = colors[0]
        elif data_prefolder == data_prefolders[1]:
            avg_data = BAPAWELD_avg_data
            std_data = BAPAWELD_std_data
            num_data = BAPAWELD_num_data
            err_data = BAPAWELD_err_data
            label = "BAPAWELD"
            color = colors[2]
        elif data_prefolder == data_prefolders[2]:
            avg_data = CBAPA_avg_data
            std_data = CBAPA_std_data 
            num_data = CBAPA_num_data 
            err_data = CBAPA_err_data 
            label = "CBAPA"
            color = colors[1]
        n_i = 0
        # print(avg_data[d_i,0,n_i,:])
        # print(avg_data[d_i,1,n_i,:])
        # ax.plot(temps,ratio_data[h_i,n_i,:],\
        ax.errorbar(avg_data[d_i,0,n_i,:],M,\
                yerr=err_data[d_i,0,n_i,:],\
                # xerr=err_data[d_i,1,n_i,:],\
                label=f"{label}",\
                color=color,\
                linestyle=styles[n_i],\
                marker='.',markersize=10,zorder=5)
        # if include_totals:
        for txt_i, txt in enumerate(M):
            ax.annotate("{:0.0f}".format(txt), (M[txt_i]+0.01, avg_data[d_i,0,n_i,txt_i]+0.01))


    # if include_totals:
    #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
    #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # ax.set_ylabel(requested_data_headers[0])
    # ax.set_xlabel(requested_data_headers[1])
    ax.set_ylabel("c/a")
    ax.set_xlabel("M")

    # ax.set_ylabel(label_from_header(header))
    # ax.set_title(f'{data_prefolder}')
    # ax.set_title('Both relaxed')
    # ax.set_xscale('log')
    # if header == requested_data_headers[-1]:
    # fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))

    #Shade the region:
    # Save current limits so autoscaling from fill doesn't move them
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()


    ax.legend()


    # plt.tight_layout()
    if save_plots:
        sav = f"{path}data/figures/c_over_a_vs_M.png"
        print(f"saving to {sav}")
        plt.savefig(sav)
    if show_plots:
        plt.show() 
    plt.close()


def gen_geometry_plot(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    attempts = [i for i in range(30)]
    N = [300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150]
    C=30
    M_sp = 20
    C_sp = [1,2,3,4,5,6,10,15,20,30,60,61]
    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    # bool_headers = [1,0,0,0,1,0,0,0,0,0]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = ["c over a"]

    data_prefolders = []

    data_prefolders.append(path + 'jobs/BAPA_')
    data_prefolders.append(path + 'jobs/BAPAWELD_')
    data_prefolders.append(path + 'jobs/CBAPA_')
    data_prefolders.append(path + 'jobs/SeqStickLognorm_')
    data_prefolders.append(path + 'jobs/DBAPA_')

    avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(M),2),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    BAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    SeqSt_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    SeqSt_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    SeqSt_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    SeqSt_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        xdata = M
        if data_prefolder == data_prefolders[4]: #DBAPA
            xdata = C_sp
        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(xdata),2),fill_value=np.nan,dtype=np.float64)
        for a_i,att in enumerate(attempts):
            n_i = 0
            n = N[0]
            size = n
            h_i = 0

            for x_i,x in enumerate(xdata):
                m = x
                m_i = x_i
                if data_prefolder == data_prefolders[2]: #CBAPA
                    c = C 
                    n = m*c
                    size = n
                if data_prefolder == data_prefolders[4]: #DBAPA
                    m = M_sp
                    c = x
                    n = m*c
                    m_i = 0

                if data_prefolder == data_prefolders[3]: #Sequential sticking
                    folder = f"{data_prefolder}{att}/N_{n}/"
                    if m != 1:
                        continue
                else:
                    folder = f"{data_prefolder}{att}/M_{m}/N_{n}/T_1000/"

                if os.path.exists(f'{folder}{n}_simData.csv'):
                    pos,radius,mass,moi = u.get_data(folder,data_index=n,relax=relax)
                    if pos is None:
                        raw_data[h_i,a_i,n_i,m_i] = np.nan
                    else:
                        effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 
                        c,b,a = u.calc_equivalent_ellipsoid_principal_axes(effective_radius,pos,mass,folder,write=False)
                        raw_data[h_i,a_i,n_i,x_i,0] = b/a
                        raw_data[h_i,a_i,n_i,x_i,1] = c/a
        
                else:
                    print(f"DNE: {folder}")


        avg_data[d_i,:,:,:,:] = np.nanmean(raw_data,axis=1)
        std_data[d_i,:,:,:,:] = np.nanstd(raw_data,axis=1)
        num_data[d_i,:,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data[d_i,:,:,:,:] = std_data[d_i,:,:,:,:]/np.sqrt(num_data[d_i,:,:,:,:])

        if data_prefolder == data_prefolders[0]:
            BAPA_avg_data = avg_data
            BAPA_std_data = std_data
            BAPA_num_data = num_data
            BAPA_err_data = err_data
        elif data_prefolder == data_prefolders[1]:
            BAPAWELD_avg_data = avg_data
            BAPAWELD_std_data = std_data
            BAPAWELD_num_data = num_data
            BAPAWELD_err_data = err_data
        elif data_prefolder == data_prefolders[2]:
            CBAPA_avg_data = avg_data
            CBAPA_std_data = std_data
            CBAPA_num_data = num_data
            CBAPA_err_data = err_data
        elif data_prefolder == data_prefolders[3]:
            SeqSt_avg_data = avg_data
            SeqSt_std_data = std_data
            SeqSt_num_data = num_data
            SeqSt_err_data = err_data
        elif data_prefolder == data_prefolders[4]:
            DBAPA_avg_data = avg_data
            DBAPA_std_data = std_data
            DBAPA_num_data = num_data
            DBAPA_err_data = err_data


    print("======================Starting figures======================")
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    # One subplot per data set
    fig, axes = plt.subplots(
        2, 2,
        figsize=(10, 10),
        sharex=True,
        sharey=True
    )

    axes = axes.flatten()

    axes_order = [0,1,2,0,3]

    for d_i, data_prefolder in enumerate(data_prefolders):

        ax = axes[axes_order[d_i]]
        X = M

        if data_prefolder == data_prefolders[0]:
            avg_data = np.copy(BAPA_avg_data)
            std_data = np.copy(BAPA_std_data)
            num_data = np.copy(BAPA_num_data)
            err_data = np.copy(BAPA_err_data) 
            label = "CAS ($N=300$)"
            color = colors[0]
            super_label = "M="
            X = M

        elif data_prefolder == data_prefolders[1]:
            avg_data = BAPAWELD_avg_data
            std_data = BAPAWELD_std_data
            num_data = BAPAWELD_num_data
            err_data = BAPAWELD_err_data
            label = "wCAS ($N=300$)"
            color = colors[2]
            super_label = "M="
            X = M

        elif data_prefolder == data_prefolders[2]:
            avg_data = CBAPA_avg_data
            std_data = CBAPA_std_data 
            num_data = CBAPA_num_data 
            err_data = CBAPA_err_data 
            label = "CNP ($C=30$)"
            color = colors[1]
            super_label = "M="
            X = M

        elif data_prefolder == data_prefolders[3]:
            avg_data = SeqSt_avg_data
            std_data = SeqSt_std_data 
            num_data = SeqSt_num_data 
            err_data = SeqSt_err_data 
            label = "Seq. Stick. ($M=1$)"
            color = colors[3]
            super_label = ""
            X = M

        elif data_prefolder == data_prefolders[4]:
            avg_data = DBAPA_avg_data
            std_data = DBAPA_std_data 
            num_data = DBAPA_num_data 
            err_data = DBAPA_err_data 
            label = "CFS ($M=20$)"
            color = colors[4]
            super_label = "C="
            X = C_sp

        n_i = 0

        ax.errorbar(
            avg_data[d_i, 0, n_i, :, 0],
            avg_data[d_i, 0, n_i, :, 1],
            yerr=err_data[d_i, 0, n_i, :, 1],
            xerr=err_data[d_i, 0, n_i, :, 0],
            label=label,
            color=color,
            # linestyle='->',
            linestyle='none',
            marker='.',
            markersize=10,
            zorder=5
        )

        for txt_i, txt in enumerate(X):
            if txt_i != np.flatnonzero(~np.isnan(avg_data[d_i, 0, n_i, :, 0]))[0]:
                final_label = ""
                offset = 0.015
                if d_i == 0:
                    offset = 0.01
            else:
                final_label = super_label
                offset = 0.01
            ax.annotate(
                "{}{:0.0f}".format(final_label,txt),
                (
                    avg_data[d_i, 0, n_i, txt_i, 0] + offset,
                    avg_data[d_i, 0, n_i, txt_i, 1] + offset
                )
            )


        #transparent BAPA points
        if axes_order[d_i] != 0:
            ax.errorbar(
                avg_data[0, 0, n_i, :, 0],
                avg_data[0, 0, n_i, :, 1],
                yerr=err_data[0, 0, n_i, :, 1],
                xerr=err_data[0, 0, n_i, :, 0],
                label="CAS ($N=300$)",
                color=colors[0],
                # linestyle='->',
                linestyle='none',
                marker='.',
                markersize=10,
                zorder=5,
                alpha=0.15
            )

        # Reference line
        ax.plot([0, 1], [0, 1], color='black', linewidth=1)

        # Shape labels
        points = [[0, 0], [1, 0], [1, 1]]
        text = ['needle', 'disk', 'sphere']
        offset = [[20, 0], [-40, 0], [-70, -20]]

        for p_i, p in enumerate(points):
            ax.plot(
                p[0], p[1],
                marker='o',
                markersize=18,
                markerfacecolor='none',
                markeredgecolor='black',
                markeredgewidth=2,
                linestyle='none',
                zorder=10,
                clip_on=False,
            )

            ax.annotate(
                text[p_i],
                xy=(p[0], p[1]),
                xytext=(offset[p_i][0], offset[p_i][1]),
                textcoords='offset points',
                ha='left',
                va='bottom',
                zorder=11
            )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        # ax.set_title(label)
        ax.legend(loc="upper left")

    # Shared axis labels
    fig.supxlabel("b/a",x=0.54)
    fig.supylabel("c/a",y=0.545)

    plt.tight_layout()

    if save_plots:
        sav = f"{path}data/figures/geometry_plot.png"
        print(f"saving to {sav}")
        plt.savefig(sav, dpi=300, bbox_inches="tight")

    if show_plots:
        plt.show()

    plt.close()



def gen_BAPA_porosity_vs_asymmetry(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    attempts = [i for i in range(30)]
    N = [300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150]
    C=30
    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    # bool_headers = [0,0,0,0,1,0,0,0,0,1,0] #Pgcs
    # bool_headers = [0,1,0,0,1,0,0,0,0,0,0] #PKBM
    bool_headers = [1,0,0,0,1,0,0,0,0,0,0] #Pabc

    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]
    requested_data_headers = list(reversed(requested_data_headers))

    data_prefolders = []

    data_prefolders.append(path + 'jobs/BAPA_')
    data_prefolders.append(path + 'jobs/CBAPA_')
    data_prefolders.append(path + 'jobs/BAPAWELD_')
    data_prefolders.append(path + 'jobs/DBAPA_')

    avg_data = np.full(shape=(len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    BAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    C_sp = [1,2,3,4,5,6,10,15,20,30,60]
    M_sp = 20
    DBAPA_avg_data = np.full(shape=(len(requested_data_headers),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    DBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    DBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    

    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        if "DBAPA" in data_prefolder:
            xdata = C_sp
        else:
            xdata = M

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            size = N[0]
            n_i = 0
            n=size

            
            for x_i,x in enumerate(xdata):

                m = x
                if "CBAPA" in data_prefolder: #CBAPA
                    size = C*m
                    n = size
                elif "DBAPA" in data_prefolder: #DBAPA
                    size = x*M_sp
                    m = M_sp
                    n = size
                                    
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_1000/"
                full_data_path = folder+f"{rel}{data_file}"
                if os.path.exists(full_data_path):
                    print(f"opening {full_data_path}")
                    with open(full_data_path,'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if size not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(size)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,n_i,x_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                else:
                    print(f"DNE: {full_data_path}")


        avg_data = np.nanmean(raw_data,axis=1)
        std_data = np.nanstd(raw_data,axis=1)
        num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
        err_data = std_data/np.sqrt(num_data)

        if data_prefolder == data_prefolders[0]:
            BAPA_avg_data = avg_data
            BAPA_std_data = std_data
            BAPA_num_data = num_data
            BAPA_err_data = err_data
        elif data_prefolder == data_prefolders[1]:
            CBAPA_avg_data = avg_data
            CBAPA_std_data = std_data
            CBAPA_num_data = num_data
            CBAPA_err_data = err_data
        elif data_prefolder == data_prefolders[2]:
            BAPAWELD_avg_data = avg_data
            BAPAWELD_std_data = std_data
            BAPAWELD_num_data = num_data
            BAPAWELD_err_data = err_data
        elif data_prefolder == data_prefolders[3]:
            DBAPA_avg_data = avg_data
            DBAPA_std_data = std_data
            DBAPA_num_data = num_data
            DBAPA_err_data = err_data
        

    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    # for t_i,temp in enumerate(temps):(NMC)

    fig,ax = plt.subplots(figsize=(8.5,7))

    for d_i,data_prefolder in enumerate(data_prefolders):
        if data_prefolder == data_prefolders[0]:
            avg_data = BAPA_avg_data
            std_data = BAPA_std_data 
            num_data = BAPA_num_data 
            err_data = BAPA_err_data 
            label = "Const. agg. size (CAS) $N=300$"
            color = colors[0]
        elif data_prefolder == data_prefolders[1]:
            avg_data = CBAPA_avg_data
            std_data = CBAPA_std_data 
            num_data = CBAPA_num_data 
            err_data = CBAPA_err_data 
            label = "Const. num. projectiles (CNP) $C=30$"
            color = colors[1]
        elif data_prefolder == data_prefolders[2]:
            avg_data = BAPAWELD_avg_data
            std_data = BAPAWELD_std_data
            num_data = BAPAWELD_num_data
            err_data = BAPAWELD_err_data
            label = "welded const. agg. size (wCAS) $N=300$"
            color = colors[2]
        elif data_prefolder == data_prefolders[3]:
            avg_data = DBAPA_avg_data
            std_data = DBAPA_std_data
            num_data = DBAPA_num_data
            err_data = DBAPA_err_data
            label = "Const. frag. size (CFS) $M=20$"
            color = colors[4]
        n_i = 0
        # print(avg_data[d_i,0,n_i,:])
        # print(avg_data[d_i,1,n_i,:])
        # ax.plot(temps,ratio_data[h_i,n_i,:],\
        ax.errorbar(avg_data[1,n_i,:],avg_data[0,n_i,:],\
                yerr=err_data[0,n_i,:],\
                xerr=err_data[1,n_i,:],\
                label=f"{label}",\
                color=color,\
                linestyle=styles[n_i],\
                marker='.',markersize=10,zorder=d_i)

        y = avg_data[0,n_i,:]
        x = avg_data[1,n_i,:]
        for i in range(len(x) - 1):
            if np.any(np.isnan([x[i], y[i], x[i+1], y[i+1]])):
                continue

            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]

            # Short arrow centered halfway between the points
            start_x = x[i]
            start_y = y[i]
            end_x   = x[i] + 0.60 * dx
            end_y   = y[i] + 0.60 * dy

            ax.annotate(
                "",
                xy=(end_x, end_y),
                xytext=(start_x, start_y),
                arrowprops=dict(
                    arrowstyle="->",
                    color=color,
                    lw=1.5
                ),
                zorder=d_i + 10
            )

        # if include_totals:
        # for txt_i, txt in enumerate(xdata):
        #     ax.annotate("{:0.0f}".format(txt), (avg_data[1,n_i,txt_i]+0.01, avg_data[0,n_i,txt_i]+0.01))



    # if include_totals:
    #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
    #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # ax.set_ylabel(requested_data_headers[0])
    # ax.set_xlabel(requested_data_headers[1])
    ax.set_ylabel(label_from_header(requested_data_headers[0]))
    ax.set_xlabel(label_from_header(requested_data_headers[1]))

    # ax.set_ylabel(label_from_header(header))
    # ax.set_title(f'{data_prefolder}')
    # ax.set_title('Both relaxed')
    # ax.set_xscale('log')
    # if header == requested_data_headers[-1]:
    # fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))

    #Shade the region:
    # Save current limits so autoscaling from fill doesn't move them
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    Pmin = 1e-12#max(xmin, 1e-12)
    Pmax = min(xmax, 1 - 1e-12)

    if Pmax <= Pmin:
        raise ValueError("No valid P range: boundary requires 0 < P < 1")

    if ax.get_xscale() == 'log':
        P = np.geomspace(Pmin, Pmax, 500)
    else:
        P = np.linspace(Pmin, Pmax, 500)

    A_boundary = 1 + 0.49 / (1 - P)**(4.0/3.0)

    ax.plot(P, A_boundary, linestyle='--', linewidth=1.5, zorder=3)

    # Fill region above the boundary up to the top of the current axes
    ax.fill_between(
        P,
        A_boundary,
        ymax,
        alpha=0.12,
        zorder=0,
        label="Allowable"
    )
    

    # Restore limits so fill doesn't change view
    ax.set_xlim(0.2, xmax)
    ax.set_ylim(ymin, ymax)

    ax.legend(loc="lower center",bbox_to_anchor=(0.5, 1.02))


    plt.tight_layout()
    if save_plots:
        sav = f"{path}data/figures/asymmetry_vs_porosity.png"
        print(f"saving to {sav}")
        plt.savefig(sav)
    if show_plots:
        plt.show() 
    plt.close()


def gen_individual_BAPA_porosity_vs_asymmetry(show_plots=True,save_plots=False,include_totals=False):
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]


    attempts = [i for i in range(30)]
    N = [300]
    M = [1,3,5,10,15,20,30,50,60,75,100,150]
    C=30
    
    data_files = []
    # data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
    data_file = "job_data.csv"


    bool_headers = [1,0,0,0,1,0,0,0,0,0]
    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

    data_prefolders = []

    data_prefolders.append(path + 'jobs/BAPA_')
    data_prefolders.append(path + 'jobs/BAPAWELD_')
    data_prefolders.append(path + 'jobs/CBAPA_')

    avg_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
    # std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    BAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    BAPAWELD_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPAWELD_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPAWELD_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # BAPAWELD_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    CBAPA_avg_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # CBAPA_std_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # CBAPA_num_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    # CBAPA_err_data = np.full(shape=avg_data.shape,fill_value=np.nan,dtype=np.float64)
    
    for d_i,data_prefolder in enumerate(data_prefolders):
        dataset_name = data_prefolder.split("/")[-1].strip("_")
        figure_folder = path+f'data/figures/{dataset_name}/'

        if save_plots and not os.path.exists(figure_folder):
            os.makedirs(figure_folder)

        relax = ("relax" in data_prefolder)
        # relax = not ("nonrelax" in data_file)
        # relax = False
        print(f"relax: {relax}")
        rel = ""
        if relax:
            rel = "relax_"

        raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(M)),fill_value=np.nan,dtype=np.float64)
        for a_i,a in enumerate(attempts):
            size = N[0]
            n_i = 0
            n=size

            
            for m_i,m in enumerate(M):
                if data_prefolder == data_prefolders[2]: #CBAPA
                    size = C*m
                    n = size
                folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_1000/"
                full_data_path = folder+f"{rel}{data_file}"
                if os.path.exists(full_data_path):
                    print(f"opening {full_data_path}")
                    with open(full_data_path,'r') as fp:
                        existing_data = fp.readlines()

                    existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
                    #even though the data can have other sizes in it, 
                    #we only want the data of size n
                    if size not in existing_sizes:
                        print(f"ERROR: Data of size {n} does not exist for {folder}.")
                        continue
                    index = existing_sizes.index(size)*4
                    existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                    existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
                    
                    for h_i,header in enumerate(requested_data_headers):
                        if header in existing_headers_for_size:
                            raw_data[h_i,a_i,n_i,m_i] = u.get_plottable_value_from_saved_value(existing_values_for_size[existing_headers_for_size.index(header)],header,folder,size,relax)
                else:
                    print(f"DNE: {full_data_path}")


        avg_data = raw_data
        # std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
        # num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
        # err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])

        if data_prefolder == data_prefolders[0]:
            BAPA_avg_data = avg_data
            # BAPA_std_data = std_data
            # BAPA_num_data = num_data
            # BAPA_err_data = err_data
        elif data_prefolder == data_prefolders[1]:
            BAPAWELD_avg_data = avg_data
            # BAPAWELD_std_data = std_data
            # BAPAWELD_num_data = num_data
            # BAPAWELD_err_data = err_data
        elif data_prefolder == data_prefolders[2]:
            CBAPA_avg_data = avg_data
            # CBAPA_std_data = std_data
            # CBAPA_num_data = num_data
            # CBAPA_err_data = err_data


    print("======================Starting figures======================")
    # print(data.shape)
    print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))

    #   plt.close("all")
    plt.rcParams.update({
        'font.size': 18,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
    })

    #Plot metric vs M for all metrics and all N and temps
    # for t_i,temp in enumerate(temps):

    fig,ax = plt.subplots()

    for d_i,data_prefolder in enumerate(data_prefolders):
        if data_prefolder == data_prefolders[0]:
            avg_data = BAPA_avg_data
            # std_data = BAPA_std_data 
            # num_data = BAPA_num_data 
            # err_data = BAPA_err_data 
            label = "BAPA"
            color = colors[0]
        elif data_prefolder == data_prefolders[1]:
            avg_data = BAPAWELD_avg_data
            # std_data = BAPAWELD_std_data
            # num_data = BAPAWELD_num_data
            # err_data = BAPAWELD_err_data
            label = "BAPAWELD"
            color = colors[2]
        elif data_prefolder == data_prefolders[2]:
            avg_data = CBAPA_avg_data
            # std_data = CBAPA_std_data 
            # num_data = CBAPA_num_data 
            # err_data = CBAPA_err_data 
            label = "CBAPA"
            color = colors[1]
        n_i = 0

        # ax.errorbar(avg_data[d_i,1,n_i,:],avg_data[d_i,0,n_i,:],\
        #         yerr=err_data[d_i,0,n_i,:],\
        #         xerr=err_data[d_i,1,n_i,:],\
        #         label=f"{label}",\
        #         color=color,\
        #         linestyle=styles[n_i],\
        #         marker='.',markersize=10,zorder=5)
        ax.plot(avg_data[1,:,n_i,:],avg_data[0,:,n_i,:],\
                label=f"{label}",\
                color=color,\
                linestyle='none',\
                marker='.',markersize=10,zorder=5)
        # if include_totals:
        # for txt_i, txt in enumerate(M):
        #     ax.annotate("{:0.0f}".format(txt), (avg_data[d_i,1,n_i,txt_i]+0.01, avg_data[d_i,0,n_i,txt_i]+0.01))


    # if include_totals:
    #   for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
    #       ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # ax.set_ylabel(requested_data_headers[0])
    # ax.set_xlabel(requested_data_headers[1])
    ax.set_ylabel(label_from_header(requested_data_headers[0]))
    ax.set_xlabel(label_from_header(requested_data_headers[1]))

    # ax.set_ylabel(label_from_header(header))
    # ax.set_title(f'{data_prefolder}')
    # ax.set_title('Both relaxed')
    # ax.set_xscale('log')
    # if header == requested_data_headers[-1]:
    # fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))

    #Shade the region:
    # Save current limits so autoscaling from fill doesn't move them
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Domain must be A>0 (avoid division by zero)
    if ax.get_xscale() == 'log':
        A = np.geomspace(max(xmin, 1e-12), xmax, 500)
    else:
        A = np.linspace(max(xmin, 1e-12), xmax, 500)

    print()


    Aminus_one = A-np.full_like(A,1)

    P_boundry = 1-np.power(0.49/Aminus_one,3.0/4.0)


    # Draw boundary
    ax.plot(A, P_boundry, linestyle='--', linewidth=1.5, zorder=3)
    # ax.plot(A, P_boundry, linestyle='--', linewidth=1.5, zorder=3, label=r"$\mathcal{P}_{abc}=1-\left(\frac{0.49}{\mathcal{A}-1}\right)^{3/4}$")

    # Shade region A > 1 + 2/P (everything above the curve)
    ax.fill_between(A, P_boundry, 0, alpha=0.12, zorder=0, label="Allowable region")

    # Restore limits so fill doesn't change view
    # ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, ymax)

    # ax.legend()


    # plt.tight_layout()
    if save_plots:
        sav = f"{path}data/figures/porosity_vs_asymmetry_BAPA.png"
        print(f"saving to {sav}")
        plt.savefig(sav)
    if show_plots:
        plt.show() 
    plt.close()



if __name__ == '__main__':
    #Do you want to see plots of the data as they are made?
    show_plots = True
    #Do you want to save the plots once they are made?
    save_plots = True
    #Do you want the number of runs next to each point on the plots
    #so you know how many more runs need to finish
    include_totals = True


    # gen_Asym_BAPA_numbers()
    # gen_agg_im_plot_BAPA(save_plots=save_plots,show_plots=show_plots)

    ##Plots for paper 2
    gen_other_BAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    gen_DBAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    gen_BAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BAPA_plots_images(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BAPA_porosity_vs_asymmetry(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_geometry_plot(save_plots=save_plots,show_plots=show_plots)
    # gen_agg_im_plot_paper2(save_plots=save_plots,show_plots=show_plots)

    # gen_BAPA_eff_rad(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_individual_BAPA_porosity_vs_asymmetry(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_c_over_a_plot(save_plots=save_plots,show_plots=show_plots)

    ##Plots for paper 1
    # gen_agg_im_plot(save_plots=save_plots,show_plots=show_plots)
    # gen_BPCA_vs_time_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_ratio_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_temp_sensitivity_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_porosity_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_gcs_csv_tables(save_plots=save_plots)
    
    


    # gen_stylized_BAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_vs_time_avg_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_seqstick_plots(distribution="lognormal")
    # gen_seqstick_plots(distribution="constant")


    # gen_relax_vs_tense_BPCA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_relax_vs_tense_seqstick_plots(distribution="lognormal",show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)

    # gen_BPCA_ratio_bugbetter_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_double_ratio_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_ratio_nonreltorel_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
    # gen_BPCA_rolling_fric_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)

    # gen_BPCA_porosity_vs_asymmetry(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)


    # T3 = 0.43736817467052647
    # T1000 = 0.3927614742047937
    # temp = 3
    # r0 = 1e-5

    # sizes=list(range(30,3001))
    # Vtot = sizes[0]*(4*np.pi/3)*(r0)**3
    # reff = (3*Vtot/(4*np.pi))**(1/3)

    # data = Tanaka(sizes,(reff/(1-T3)**(1/3))/(5/3)**(1/2),temp)
    
    # fig,ax = plt.subplots()
    # ax.plot(sizes,data,temp,\
    #       label=f"Tanaka prediction T={temp}")

    # ax.set_ylim(0.42,0.46)

    
    # plt.show() 
    # plt.close()