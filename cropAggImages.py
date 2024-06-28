import glob as g
from PIL import Image
import os
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	data_directory = input_json["data_directory"]
	job_group = "const_rel"
	job_group = "lognorm_rel"

	image_path = data_directory + "data/figures/aggRenders/"
	save_path = image_path+"edited/"
	
	Nums = [30,100]
	Nums = [300]
	temps = [3,10,30,100,300,1000]

	for N in Nums:
		for t in temps:
			# N=300
			# t=3
			glob_me = image_path+f'agg-{job_group}*_a-*_N-{N}_T-{t}.png'
			print(glob_me)
			image = g.glob(glob_me)[0]
			
			im = Image.open(image)
			width,height = im.size

			left = width*1/5
			top = 0
			right = width*4/5
			bottom = height
			

			im1 = im.crop((left,top,right,bottom))
			# im1.show()
			# exit(0)
			whole_path = image.split('/')
			crop_name = whole_path[-1][:-4]+'_cropped'+whole_path[-1][-4:]
			# image = '/'.join(path)

			im1.save(save_path+crop_name)

if __name__ == '__main__':
	main()