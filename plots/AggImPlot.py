import matplotlib.pyplot as plt
import glob as g
import numpy as np
import os
import json
from PIL import Image
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.ticker as plticker
from matplotlib.patches import FancyArrowPatch


relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

plt.rcdefaults()
plt.rcParams['font.size'] = 30

def main():
	constant = 1.5
	fig,ax = plt.subplots(figsize=(constant*1152*6/500 - .10, constant*1080*3/500 ))

	grid = ImageGrid(fig,111,nrows_ncols=(3,6),axes_pad=0,share_all=True)

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	job_group = "const_rel"
	job_group = "lognormrelax"

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
			print(glob_me)
			image = g.glob(glob_me)[0]
			im = Image.open(image)

			# Convert the image to a NumPy array
			image_array = np.array(im)

			# Print the shape of the array to understand its structure
			# print(f"Array shape: {image_array.shape}")

			# Check if the image is grayscale
			if image_array.ndim == 2:
				cmaps.append('gray')
			else:
				cmaps.append(None)
				# print("The image is grayscale, setting cmap to 'grey'.")
			# elif image_array.ndim == 3 and (image_array[:,:,0] == image_array[:,:,1]).all() and (image_array[:,:,1] == image_array[:,:,2]).all():
				# print("The image is RGB but all channels are identical, which means it's effectively grayscale.")

			# print(im.size)
			images.append(im)

	# print(images)
	# grid[0].imshow(images[0])
	# plt.xticks(np.arange(0,1,step=1/6))
	ax.xaxis.set(ticks=[5,15,25,35,45,55,60],
				ticklabels=['3','10','30','100','300','1000',''],
				label="x")
	# loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
	# ax.xaxis.set_major_locator(loc)
	ax.yaxis.set(ticks=[5,15,25,30],ticklabels=['300','100','30',''])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	index = 0
	for axe,im in zip(grid,images):
		axe.axis('off')
		axe.imshow(im, cmap=cmaps[index])
		index

		# axe.imshow(im,cmap=None)
	# axe.imshow(im, cmap=None)
	ax.set_xlabel('Temp (K)')
	ax.set_ylabel('Number of Particles')


	h = 0.61
	start = (0.35, h)    # figure-fraction coords
	end   = (0.68, h)

	arrow = FancyArrowPatch(
	    start, end,
	    transform=fig.transFigure,
	    arrowstyle='-|>',
	    mutation_scale=25,
	    lw=5, color='black',
	    zorder=200
	)
	fig.add_artist(arrow)

	fig.text(
	    0.515, end[1]+0.03,       # slightly above the head
	    "Density increases with temperature",
	    transform=fig.transFigure,
	    color='black',
	    fontsize=30,
	    ha='center'
	)


	# plt.savefig(path + f'data/figures/ColoredAggComp_{job_group}.png',dpi=1500)
	plt.show()

if __name__ == '__main__':
	main()