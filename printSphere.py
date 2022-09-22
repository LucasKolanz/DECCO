import numpy as np
import matplotlib.pyplot as plt


def main():
	file = "enclosingSphere.txt"
	with open(file, 'r') as fp:
		lines = fp.readlines()

	lines = [i.strip('\n').split(',') for i in lines]
	# lines[:][2] = [i.strip('\n') for i in lines[:][2]]
	points = np.float_(lines[:-2])
	center = lines[-2]
	center[0] = center[0].strip("center: ")
	center = [float(i) for i in center]
	print(center)
	radius = float(lines[-1][0].strip('radius: ')) 
	print(radius)

	# circlexy = plt.Circle((center[0],center[1]),radius, fill=False)
	# circleyz = plt.Circle((center[1],center[2]),radius, fill=False)
	# circlexz = plt.Circle((center[0],center[2]),radius, fill=False)

	# fig,ax = plt.subplots(2,2)
	# #xy plane
	# ax[0,0].scatter(points[:,0],points[:,1],color='b')
	# ax[0,0].add_patch(circlexy)
	# #yz plane
	# ax[0,1].scatter(points[:,1],points[:,2],color='b')
	# ax[0,1].add_patch(circleyz)
	# #xz plane
	# ax[1,0].scatter(points[:,0],points[:,2],color='b')
	# ax[1,0].add_patch(circlexz)

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')

	# ax.scatter(points[:,0],points[:,1],points[:,2])
	init_rad = 1e-5
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	for pt in points:
		print(pt)
		x = init_rad*(np.cos(u)*np.sin(v)) + pt[0]
		y = init_rad*(np.sin(u)*np.sin(v)) + pt[1]
		z = init_rad*(np.cos(v)) + pt[2]
		ax.plot_wireframe(x, y, z, color="b")

	#u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = radius*(np.cos(u)*np.sin(v)) + center[0]
	y = radius*(np.sin(u)*np.sin(v)) + center[1]
	z = radius*(np.cos(v)) + center[2]
	ax.plot_wireframe(x, y, z, color="r")

	plt.show()


if __name__ == '__main__':
	main()