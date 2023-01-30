import glob as g
from PIL import Image

def main():
	image_path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/agg_pics/'
	Nums = [30,100,300]
	temps = [3,10,30,100,300,1000]

	for N in Nums:
		for t in temps:
			# N=300
			# t=3
			image = g.glob(image_path+'originals/'+'N{}T{}A*_1.png'.format(N,t))[0]
			
			im = Image.open(image)
			width,height = im.size

			left = width*1/5
			top = 0
			right = width*4/5
			bottom = height
			

			im1 = im.crop((left,top,right,bottom))
			# im1.show()
			# exit(0)
			path = image.split('/')
			path[-2] = 'edited'
			path[-1] = path[-1][:-4]+'cropped'+path[-1][-4:]
			image = '/'.join(path)
			# print(image)
			# exit(0)
			im1.save(image)

if __name__ == '__main__':
	main()