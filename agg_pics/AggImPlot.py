import matplotlib.pyplot as plt
import glob as g
import numpy as np
from PIL import Image
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.ticker as plticker


def main():
	fig,ax = plt.subplots(figsize=(479*6/200,368*3/200))
	grid = ImageGrid(fig,111,nrows_ncols=(3,6),axes_pad=0,share_all=True)
	
	image_path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/agg_pics/'
	Nums = [30,100,300]
	temps = [3,10,30,100,300,1000]

	images = []
	for N in Nums:
		for t in temps:
			# image = g.glob(image_path+'edited/N{}T{}A*'.format(N,t))[0]
			image = g.glob(image_path+'edited/N{}T{}A*_1*.png'.format(N,t))[0]
			im = Image.open(image)
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
	for axe,im in zip(grid,images):
		axe.axis('off')
		axe.imshow(im)

	ax.set_xlabel('Temp (K)')
	ax.set_ylabel('Aggregate Size')
	plt.savefig('../figures/AggComp.png')
	plt.show()

if __name__ == '__main__':
	main()