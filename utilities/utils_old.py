import numpy as np
import random
from scipy.spatial.transform import Rotation as R
import os,glob
import sys
import json
import matplotlib.pyplot as plt
#from treelib import Node, Tree
import time
from itertools import combinations
# cwd = os.getcwd()
# os.system("cd /home/kolanzl/Open3D/build")
# sys.path.append("/home/kolanzl/Open3D/build/lib/python_package/open3d")
import open3d as o3d
##include <pybind11/stl.h>`? Or <pybind11/complex.h>,
# <pybind11/functional.h>, <pybind11/chrono.h>

data_columns = 11

#next three functions translated from matlab code from
#https://blogs.mathworks.com/cleve/files/menger.m
# def menger(level):
# 	V = [[-3,-3,-3],[-3,-3,3],[-3,3,-3],[-3,3,3],[3,-3,-3],[3,-3,3],[3,3,-3],[3,3,3]]
# 	V = np.array(V) 
# 	V = sponge(V,level)
# 	return V

# def sponge(V,level):

# 	if level > 0:
# 		V = V/3
# 		for x in [-2,0,2]:
# 			for y in [-2,0,2]:
# 				for z in [-2,0,2]:
# 					if np.sum(np.array([x,y,z])) > 0:
# 						sponge(V)
# 	# else:
# 		# cube(V)			
# 	return V	

# def cube(V):
	# return



def plot(verts,center,radius):
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	# print('')
	# print(verts)
	# print(verts[:][0])
	# exit(0)
	for i in range(len(verts)):
		print(verts[i],i)
		ax.scatter(verts[i][0],verts[i][1],verts[i][2],marker='*',color='b')
	ax.scatter(center[0],center[1],center[2],marker='.',color='r')
	ax.set_xlabel('X (row)')
	ax.set_ylabel('Y (col)')
	ax.set_zlabel('Z (dep)')

	# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	# x = np.cos(u)*np.sin(v) - center[0]
	# y = np.sin(u)*np.sin(v) - center[1]
	# z = np.cos(v) - center[2]
	# ax.plot_wireframe(x, y, z, color="r")
	plt.show()

def get_data_file(data_folder,data_index=-1):
	files = os.listdir(data_folder)

	try:
		file_indicies = np.array([file.split('_')[0] for file in files\
					if file.endswith("simData.csv")],dtype=np.int64)
 
	except: 
		files = [file for file in files if file.endswith('simData.csv')]
		files = [file for file in files if '_' in file]
		file_indicies = np.array([int(file.split('_')[0]) for file in files],dtype=np.int64)
	# 	file_indicies = 

	if data_index == -1:
		index = np.max(file_indicies)
	else:
		index = data_index

	# print("index: {}".format(index))

	data_file = [file for file in files \
				if file.endswith("simData.csv") and file.startswith(str(index))]

	if len(data_file) == 1:
		return data_file[0]
	elif len(data_file) == 2 and len(set(file_indicies)) == 1:
		if (len(data_file[0]) > len(data_file[1])):
			return data_file[0] 
		else:
			return data_file[1]
	else:
		data_file = [file for file in files \
				if file.endswith("simData.csv") and file.startswith(str(index)+'_2')]
		# print(files)
		if len(data_file) == 1:
			return data_file[0]
		elif len(data_file) == 2:
			if len(data_file[0]) > len(data_file[1]):
				return data_file[0]
			else:
				return data_file[1]
		print("data file in folder '{}' not found.".format(data_folder))
		print("Now exiting.")
		exit(-1)

def get_energy_file(data_folder,data_index=-1):
	files = os.listdir(data_folder)
	# print(files)
	try:
		file_indicies = np.array([file.split('_')[0] for file in files\
					if file.endswith("energy.csv")],dtype=np.int64)
	except: 
		files = [file for file in files if file.endswith('energy.csv')]
		files = [file for file in files if '_' in file]
		file_indicies = np.array([int(file.split('_')[0]) for file in files],dtype=np.int64)
	# 	file_indicies = 

	if data_index == -1:
		index = np.max(file_indicies)
	else:
		index = data_index

	# print(file_indicies)
	# print(np.max(file_indicies))

	# print("index: {}".format(index))

	data_file = [file for file in files \
				if file.endswith("energy.csv") and file.startswith(str(index))]
	# print(data_file)

	if len(data_file) == 1 or len(set(data_file)) == 1:
		return data_file[0]
	elif len(data_file) == 2 and len(set(file_indicies)) == 1:
		if (len(data_file[0]) > len(data_file[1])):
			return data_file[0] 
		else:
			return data_file[1]
	else:
		data_file = [file for file in files \
				if file.endswith("energy.csv") and file.startswith(str(index)+'_2')]
		if len(data_file) == 1:
			return data_file[0]
		elif len(data_file) == 2:
			if len(data_file[0]) > len(data_file[1]):
				return data_file[0]
			else:
				return data_file[1]
		print("energy file in folder '{}' not found.".format(data_folder))
		print("Now exiting.")
		exit(-1)

def get_last_line_data(data_folder,data_index=-1):
	# data_headers = np.loadtxt(data_folder + data_file,skiprows=0,dtype=str,delimiter=',')[0]
	data_file = get_data_file(data_folder,data_index)
	# print("data file: {}".format(data_file))
	try:
		data = np.loadtxt(data_folder + data_file,skiprows=1,dtype=float,delimiter=',')
		if data.ndim > 1:
			data = data[-1]
		# print(data)
		# print(data_folder + data_file)
	except Exception as e:
		with open(data_folder + data_file) as f:
			for line in f:
				pass
			last_line = line
		data = np.array([last_line.split(',')],dtype=np.float64)
		# print(data)
		print("ERROR CAUGHT getting data in folder: {}".format(data_folder))
		print(e)
		return None
	# print("DATA LEN: {} for file {}{}".format(data.size,data_folder,data_file))
	# print("FOR {} Balls".format(data.size/11))
	return format_data(data)

def get_last_line_energy(data_folder,data_index=-1):
	energy_file = get_energy_file(data_folder,data_index)
	try:
		energy = np.loadtxt(data_folder + energy_file,skiprows=1,dtype=float,delimiter=',')
		if energy.ndim > 1:
			energy = energy[-1]
		# print(energy)
	except Exception as e:
		with open(data_folder + energy_file) as f:
			for line in f:
				pass
			last_line = line
		energy = np.array([last_line.split(',')],dtype=np.float64)
		print("ERROR CAUGHT getting energy in folder: {}".format(data_folder))
		print(e)
		# return None

	# print("DATA LEN: {} for file {}{}".format(data.size,data_folder,data_file))
	# print("FOR {} Balls".format(data.size/11))
	return energy

def get_constants(data_folder,data_index=-1):
	data_file = get_data_file(data_folder,data_index)
	data_file = data_file.replace('simData','constants')	
	data_constants = np.loadtxt(data_folder+data_file,skiprows=0,dtype=float,delimiter=',')[0]
	return data_constants[0],data_constants[1],data_constants[2]

def get_all_constants(data_folder,data_index=-1):
	data_file = get_data_file(data_folder,data_index)
	data_file = data_file.replace('simData','constants')	
	data_constants = np.loadtxt(data_folder+data_file,skiprows=0,dtype=float,delimiter=',')
	return data_constants

def format_data(data):
	data = np.reshape(data,(int(data.size/data_columns),data_columns))
	data = data[:,:3]
	return data

def COM(data_folder,data_index=-1):
	data = get_last_line_data(data_folder,data_index)
	consts = get_all_constants(data_folder,data_index)
	com = np.array([0,0,0],dtype=np.float64)
	mtot = 0

	for ball in range(data.shape[0]):
		com += consts[ball][0]*data[ball]
		mtot += consts[ball][0]

	return mtot*com

def get_data(data_folder,data_index=-1):
	if data_folder == '/home/kolanzl/Desktop/bin/merger.csv':
		data = np.loadtxt(data_folder,delimiter=',')
		radius = 1
		mass = 1
		moi = 1
	else:
		data_file = get_data_file(data_folder,data_index)
		radius,mass,moi = get_constants(data_folder,data_index)
		data = get_last_line_data(data_folder,data_index)
	return data,radius,mass,moi

def get_data_range(data_folder,data_index=-1):
	if data_folder == '/home/kolanzl/Desktop/bin/merger.csv':
		data = np.loadtxt(data_folder,delimiter=',')
		radius = 1
		mass = 1
		moi = 1
	else:
		data = get_last_line_data(data_folder,data_index)
		radius,m,moi = get_constants(data_folder,data_index)

	max_x = np.max(data[:,0]) + radius
	min_x = np.min(data[:,0]) - radius
	max_y = np.max(data[:,1]) + radius
	min_y = np.min(data[:,1]) - radius
	max_z = np.max(data[:,2]) + radius
	min_z = np.min(data[:,2]) - radius
	
	return max_x,min_x,max_y,min_y,max_z,min_z

#following functions taken from 
#http://www.open3d.org/docs/release/tutorial/geometry/voxelization.html#Voxel-carving
def xyz_spherical(xyz):
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]
	r = np.sqrt(x * x + y * y + z * z)
	r_x = np.arccos(y / r)
	r_y = np.arctan2(z, x)
	return [r, r_x, r_y]


def get_rotation_matrix(r_x, r_y):
	rot_x = np.asarray([[1, 0, 0], [0, np.cos(r_x), -np.sin(r_x)],
						[0, np.sin(r_x), np.cos(r_x)]])
	rot_y = np.asarray([[np.cos(r_y), 0, np.sin(r_y)], [0, 1, 0],
						[-np.sin(r_y), 0, np.cos(r_y)]])
	return rot_y.dot(rot_x)


def get_extrinsic(xyz):
	rvec = xyz_spherical(xyz)
	r = get_rotation_matrix(rvec[1], rvec[2])
	t = np.asarray([0, 0, 2]).transpose()
	trans = np.eye(4)
	trans[:3, :3] = r
	trans[:3, 3] = t
	return trans


def preprocess(model):
	min_bound = model.get_min_bound()
	max_bound = model.get_max_bound()
	center = min_bound + (max_bound - min_bound) / 2.0
	scale = np.linalg.norm(max_bound - min_bound) / 2.0
	# scale = 1
	vertices = np.asarray(model.vertices)
	vertices -= center
	model.vertices = o3d.utility.Vector3dVector(vertices / scale)
	return model

def preprocess_pt(model):
	min_bound = model.get_min_bound()
	max_bound = model.get_max_bound()
	print("max: {}\tmin: {}".format(max_bound,min_bound))
	center = min_bound + (max_bound - min_bound) / 2.0
	scale = np.linalg.norm(max_bound - min_bound) / 2.0
	# scale = 1
	vertices = np.asarray(model.points)
	vertices -= center
	model.points = o3d.utility.Vector3dVector(vertices / scale)
	return model


def vox_carve(mesh,
				  cubic_size,
				  voxel_resolution,
				  w=300,
				  h=300,
				  use_depth=True,
				  surface_method='pointcloud'):
	# mesh.compute_vertex_normals()
	mesh.estimate_normals()
	camera_sphere = o3d.geometry.TriangleMesh.create_sphere()

	# setup dense voxel grid
	voxel_carving = o3d.geometry.VoxelGrid.create_dense(
		width=cubic_size,
		height=cubic_size,
		depth=cubic_size,
		voxel_size=cubic_size / voxel_resolution,
		origin=[-cubic_size / 2.0, -cubic_size / 2.0, -cubic_size / 2.0],
		color=[1.0, 0.7, 0.0])
	print("Vox size: {}".format(cubic_size / voxel_resolution))
	# rescale geometry
	camera_sphere = preprocess(camera_sphere)
	mesh = preprocess_pt(mesh)

	# setup visualizer to render depthmaps
	vis = o3d.visualization.Visualizer()
	vis.create_window(width=w, height=h, visible=False)
	vis.add_geometry(mesh)
	vis.get_render_option().mesh_show_back_face = True
	ctr = vis.get_view_control()
	param = ctr.convert_to_pinhole_camera_parameters()

	# carve voxel grid
	pcd_agg = o3d.geometry.PointCloud()
	centers_pts = np.zeros((len(camera_sphere.vertices), 3))
	for cid, xyz in enumerate(camera_sphere.vertices):
		# get new camera pose
		trans = get_extrinsic(xyz)
		param.extrinsic = trans
		c = np.linalg.inv(trans).dot(np.asarray([0, 0, 0, 1]).transpose())
		centers_pts[cid, :] = c[:3]
		ctr.convert_from_pinhole_camera_parameters(param)

		# capture depth image and make a point cloud
		vis.poll_events()
		vis.update_renderer()
		depth = vis.capture_depth_float_buffer(False)
		pcd_agg += o3d.geometry.PointCloud.create_from_depth_image(
			o3d.geometry.Image(depth),
			param.intrinsic,
			param.extrinsic,
			depth_scale=1)

		# depth map carving method
		if use_depth:
			voxel_carving.carve_depth_map(o3d.geometry.Image(depth), param)
		else:
			voxel_carving.carve_silhouette(o3d.geometry.Image(depth), param)
		# print("Carve view %03d/%03d" % (cid + 1, len(camera_sphere.vertices)))
	vis.destroy_window()

	# add voxel grid survace
	print('Surface voxel grid from %s' % surface_method)
	if surface_method == 'pointcloud':
		voxel_surface = o3d.geometry.VoxelGrid.create_from_point_cloud_within_bounds(
			pcd_agg,
			voxel_size=cubic_size / voxel_resolution,
			min_bound=(-cubic_size / 2, -cubic_size / 2, -cubic_size / 2),
			max_bound=(cubic_size / 2, cubic_size / 2, cubic_size / 2))
	elif surface_method == 'mesh':
		voxel_surface = o3d.geometry.VoxelGrid.create_from_triangle_mesh_within_bounds(
			mesh,
			voxel_size=cubic_size / voxel_resolution,
			min_bound=(-cubic_size / 2, -cubic_size / 2, -cubic_size / 2),
			max_bound=(cubic_size / 2, cubic_size / 2, cubic_size / 2))
	else:
		raise Exception('invalid surface method')
	voxel_carving_surface = voxel_surface + voxel_carving

	return voxel_carving_surface, voxel_carving, voxel_surface

def dist(pt1,pt2):
	return np.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2)


class datamgr(object):
	"""docstring for datamgr"""
	# def __init__(self, data_folder,voxel_buffer=5,ppb=3000):
	def __init__(self, data_folder,index=-1,ppb=30000,Temp=-1):
		super(datamgr, self).__init__()
		self.data_folder = data_folder
		self.index = index
		if data_folder != '/home/kolanzl/Desktop/bin/merger.csv' and Temp < 0:
			self.Temp = int(self.data_folder.split('/')[-2].strip('T_'))
		else:
			self.Temp = Temp
		#how many points in single ball pointcloud shell
		self.ppb = ppb
		self.data,self.radius,self.mass,self.moi = get_data(self.data_folder,self.index)
		self.nBalls = self.data.shape[0]
		# self.buffer = voxel_buffer # how many extra voxels in each direction 
		self.data_range = get_data_range(self.data_folder,self.index)

	def shift_to_first_quad(self,data_range=None):
		if data_range is None:
			data_range = get_data_range(self.data_folder,self.index)
		# print("SHIFTED")

		self.data[:,0] -= data_range[1] 
		self.data[:,1] -= data_range[3] 
		self.data[:,2] -= data_range[5] 


	def vox_init(self,num_vox):
		data_abs_max = max(self.data_range,key=abs) 
		self.vox_size = (data_abs_max*2)/num_vox
		# self.vox_size = 
		# print(self.vox_per_radius)
		# self.vox_rep = np.zeros((num_vox+self.buffer*2,num_vox+self.buffer*2,num_vox+self.buffer*2))

	#Function written by chatGPT
	def rotation_matrix(v1, v2):
		"""
		Returns the rotation matrix between two vectors v1 and v2.
		Both v1 and v2 must be numpy arrays with the same shape.

		:param v1: First vector
		:param v2: Second vector
		:return: Rotation matrix
		"""
		v1 = np.array(v1)
		v2 = np.array(v2)
		if v1.shape != v2.shape:
			raise ValueError("Both vectors must have the same shape.")
		v1 = v1 / np.linalg.norm(v1)
		v2 = v2 / np.linalg.norm(v2)
		v = np.cross(v1, v2)
		s = np.linalg.norm(v)
		c = np.dot(v1, v2)
		vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
		rotation_matrix = np.eye(v1.shape[0]) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
		return rotation_matrix

	def orient_data(self):
		max_lengsq = -1
		pt1 = []
		pt2 = []
		for i,p1 in enumerate(self.data):
			for j,p2 in enumerate(self.data):
				if i != j:
					lengsq = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
					if max_lengsq < lengsq:
						max_lengsq = lengsq
						pt1 = p1
						pt2 = p2
		print(max_lengsq)
		print(pt1)
		print(pt2)

	def gen_whole_pt_cloud(self):
		# self.orient_data()
		# exit(0)
		self.shift_to_first_quad()

		radii = np.linspace(self.radius/100,self.radius,100)

		accum = [self.ppb*(radius**2/self.radius**2) for radius in radii]
		accum = np.array(accum,dtype=int)
		accum = np.where(accum < 100, 100, accum)

		return_array = np.zeros((self.data.shape[0]*np.sum(np.array(accum)),3))
		for ind,pt in enumerate(self.data):
			for i,radius in enumerate(radii):
				start_index = int(ind*np.sum(accum)) + np.sum(accum[:i])
				end_index = int(ind*np.sum(accum)) + np.sum(accum[:i+1])
				return_array[start_index:end_index] = self.gen_pt_cloud(pt,radius,accum[i])
		return return_array

	#evenly spaced points on sphere code form:
	#http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
	def gen_pt_cloud(self,pt,radius,num_pts):
		goldenRatio = (1 + 5**0.5)/2

		return_array = np.zeros((num_pts,3))	
		i = np.arange(0, num_pts)

		theta = 2 * np.pi * i / goldenRatio
		phi = np.arccos(1 - 2*(i+0.5)/num_pts)

		# return_array = np.zeros((self.ppb,3))
		return_array[:,0] = radius*(np.cos(theta) * np.sin(phi)) + pt[0]
		return_array[:,1] = radius*(np.sin(theta) * np.sin(phi)) + pt[1]
		return_array[:,2] = radius*np.cos(phi) + pt[2]
		
		return return_array

	def get_center_pt(self,ind):
		return [ind[0]*self.radius+self.radius/2,ind[1]*self.radius+self.radius/2,ind[2]*self.radius+self.radius/2]


class o3doctree(object):
	"""docstring for o3doctree"""
	def __init__(self, data_folder=None,ppb=30000,verbose=False,overwrite_data=False, \
				visualize_pcd=False, visualize_octree=False, index=-1,Temp=-1):
	# def __init__(self, data_folder, max_depth=8,ppb=600000,verbose=False):
		super(o3doctree, self).__init__()
		self.data_folder = data_folder
		self.ppb = ppb
		self.verbose = verbose
		self.overwrite_data = overwrite_data
		self.visualize_pcd = visualize_pcd
		self.visualize_octree = visualize_octree
		self.bestfitlen = 4
		self.index = index
		if Temp > 0:
			self.Temp = Temp


	def make_tree(self):

		self.dm = datamgr(self.data_folder,self.index,self.ppb,Temp=self.Temp)

		bounds = [self.dm.data_range[0]-self.dm.data_range[1],self.dm.data_range[2]-self.dm.data_range[3],self.dm.data_range[4]-self.dm.data_range[5]]
		max_bound = max(bounds) + 2*self.dm.radius
		n=0
		rad = 999
		while rad > self.dm.radius:
			n += 1
			rad = max_bound/(2**n)

		self.max_depth = n


		pcd_file = self.data_folder + "pointcloud_ppb-{}.pcd".format(self.ppb)
		oct_file = self.data_folder + "octree_ppb-{}".format(self.ppb)
		fractdim_data_file = self.data_folder + "fractdim_ppb-{}.csv".format(self.ppb)
		
		make_data = False
		try:
			print("Loading FD data for :{}".format(fractdim_data_file))
			d = np.loadtxt(fractdim_data_file, delimiter=',')
			# print(d)
			# print(np.sum(d[:,1]))
			if np.sum(d[:,1]) == 0:
				print("Data needs recomputing")
				os.remove(fractdim_data_file)
				make_data = True
			else:
				with open(fractdim_data_file,'r') as f:
					header = f.readline()
				self.octree_size = float(header.strip('\n').strip('# '))
				self.s_data = d[:,0]
				self.Ns_data = d[:,1]

			# with open(fractdim_data_file,'r') as f:
			# 	data = f.readlines()
			# 	print(data)
			# if len(data) == 0:
			# 	make_data = True
			# if np.sum(d[:,0]) == 0:
		# except IOError:
		except:
			# print("Computing data for :{}".format(fractdim_data_file))
			make_data = True

		octree = []
		if make_data or self.overwrite_data:
			octverb = ''
			print("Generating FD data for :{}".format(fractdim_data_file))
			if os.path.isfile(oct_file) and not self.overwrite_data:
				octstart = time.process_time()
				octverb = 'Getting'
				octree = o3d.io.read_octree(oct_file)
				if self.visualize_octree:
					self.show_octree(self.verbose)
			else:
				if self.verbose:
					pcdstart = time.process_time()
				
				pcdverb = ''
				# if os.path.isfile(pcd_file):
				pcd = []
				make_pcd_data = False
				try:
					with open(pcd_file,'r') as f:
						data = f.readlines()
						if len(data) == 0:
							make_pcd_data = True
							os.remove(pcd_file)
				# except IOError:
				except:
					make_pcd_data = True
				if not make_pcd_data and not self.overwrite_data:
					pcdverb = 'Getting'
					pcd = o3d.io.read_point_cloud(pcd_file)
					if self.visualize_pcd:
						self.show_pcd(pcd,self.verbose)
				else:
					pcdverb = 'Making'
					# radii = np.linspace(self.dm.radius/100,self.dm.radius,100)
					# accum = [self.dm.ppb*(radius**2/self.dm.radius**2) for radius in radii]
					# accum = np.array(accum,dtype=int)
					# accum = np.where(accum < 100, 100, accum)
					point_cloud = self.dm.gen_whole_pt_cloud()
					# point_cloud = point_cloud[:sum(accum)]
					pcd = o3d.geometry.PointCloud()
					pcd.points = o3d.utility.Vector3dVector(point_cloud)
					pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=point_cloud.shape))
					o3d.io.write_point_cloud(pcd_file, pcd)
					if self.visualize_pcd:
						self.show_pcd(pcd,self.verbose)
					# exit(0)
			
				if self.verbose:
					pcdend = time.process_time()
					print("{} pcd took {:.2f}s".format(pcdverb,pcdend-pcdstart))

				octverb = 'Making'
				octstart = time.process_time()

				octree = o3d.geometry.Octree(max_depth=self.max_depth)
				if self.visualize_octree:
					self.show_octree(octree,self.verbose)
				octree.convert_from_point_cloud(pcd, size_expand=0.01)

				o3d.io.write_octree(oct_file, octree)


			if self.verbose:
				octend = time.process_time()
				print("{} octree took {:.2f}s".format(octverb, octend-octstart))
				start = time.process_time()
				print("Starting octree traversal")
		
			self.tree_info = []
			for i in range(self.max_depth):
				self.tree_info.append(0)
			# print(self.tree_info)
			if self.verbose:
				end = time.process_time()
				print("Traversing octree took {:.2f}s".format(end-start))
			octree.traverse(self.f_traverse)

			self.s_data = np.zeros((self.max_depth))
			self.Ns_data = np.zeros((self.max_depth))
			for i in range(self.max_depth):
				self.s_data[i] = (2**-(i+1))
				self.Ns_data[i] = self.tree_info[i]
			save_data = np.zeros((self.max_depth,2))
			save_data[:,0] = self.s_data
			save_data[:,1] = self.Ns_data
			np.savetxt(fractdim_data_file,save_data, delimiter=',',header=str(octree.size))
			self.octree_size = octree.size
		# else:
			
	#TODO  This function should find the orientation that minimizes 
	#	   the original fractal dimension (depth of 1)
	# def point_orientation(self,point_cloud):
	# 	# print(point_cloud)
	# 	# exit(0)
	# 	i = 0
	# 	best_i = 0
	# 	rotations = []
	# 	pcd = o3d.geometry.PointCloud()
	# 	while (i < 10):
	# 		xrot = random.uniform(0,360)
	# 		yrot = random.uniform(0,360)
	# 		zrot = random.uniform(0,360)
	# 		rotation_matrix = R.from_euler('xyz',[xrot,yrot,zrot],degrees=True).as_matrix()
	# 		rotations.append(rotation_matrix)
	# 		temp_point_cloud[:] = rotation_matrix @ point_cloud[:]
	# 		# temp_point_cloud = [rotation_matrix @ i for i in point_cloud]
	# 		pcd.points = o3d.utility.Vector3dVector(temp_point_cloud)
	# 		octree = o3d.geometry.Octree(max_depth=1)#check max_depth def
	# 		octree.convert_from_point_cloud(pcd, size_expand=0.01)#check size_expand def
	# 		self.tree_info = [0]
			
	# 		octree.traverse(self.f_traverse)
	# 		print(self.tree_info)
	# 		exit(0)
	# 		i+=1
	# 		# print(rotation_matrix.as_matrix())
	# 	exit(0)
	# 	return o3d.utility.Vector3dVector(point_cloud)

	# def add_menger_points(self,data):
	# 	dlen = data.shape
	# 	print(dlen)
	# 	exit(0)

	def test_menger_sponge(self):
		self.data_folder = '.'
		merger_file = '/home/kolanzl/Desktop/bin/merger.csv'
		self.dm = datamgr(merger_file,self.ppb)

		self.tree_info = []
		max_depth = 8
		self.s_data = np.zeros((max_depth))
		self.Ns_data = np.zeros((max_depth))
		for i in range(max_depth):
			self.tree_info.append(0)
		menger_points = np.loadtxt(merger_file,delimiter=',')
		# menger_points = self.add_menger_points(menger_points)

		pcd = o3d.geometry.PointCloud()
		pcd.points = o3d.utility.Vector3dVector(menger_points)
		pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=menger_points.shape))

		octree = o3d.geometry.Octree(max_depth=max_depth)
		octree.convert_from_point_cloud(pcd, size_expand=0.01)
		# o3d.visualization.draw_geometries([octree])
		# exit(0)

		octree.traverse(self.f_traverse)
		for i in range(max_depth):
			self.s_data[i] = (2**-(i+1))
			self.Ns_data[i] = self.tree_info[i]
		
		self.calc_fractal_dimension(True)
		# self.Ns_data[-1] -=
		print(np.log(self.Ns_data)/np.log(1/self.s_data))


	def show_octree(self,octree,verbose):
		if verbose:
			start = time.process_time()
		o3d.visualization.draw_geometries([octree])
		if verbose:
			end = time.process_time()
			print("Visualizing octree took {}".format(end-start))

	def show_pcd(self,pcd,verbose):
		if verbose:
			start = time.process_time()
		o3d.visualization.draw_geometries([pcd])
		if verbose:
			end = time.process_time()
			print("Visualizing pcd took {}".format(end-start))	

	def bestfit(self,x_data,y_data,length,min_rang=None,max_range=None):
		acceptable_indicies = []
		if min_rang is not None and max_range is not None:
			for i in range(len(x_data)):
				#Note that min_range/max_range is before taking the inverse (1/min_range etc)
				#so min_range is actually larger than max_range. 
				if x_data[i] < min_rang and x_data[i] > max_range:
					acceptable_indicies.append(i)

			x_data = x_data[acceptable_indicies]
			y_data = y_data[acceptable_indicies]
		if length > len(x_data):
			length = len(x_data)
		if length == 1:
			print("ERROR: Cannot fit to one data point")

		index_list = np.arange(0,len(x_data))
		rsq = 0
		rsq_index = 0 
		fit = []
		combos = combinations(index_list,length)
		combos = [list(c) for c in list(combos)]
		winning_combo = []

		for i,comb in enumerate(combos):
			# print(comb)
			ind_combo = np.array(sorted(comb))
			fit_x_pts = np.log(1/x_data[ind_combo])
			fit_y_pts = np.log(y_data[ind_combo])
			fit.append(np.polyfit(fit_x_pts,fit_y_pts,1))
			y_predict = np.array(fit_x_pts*fit[-1][0] + fit[-1][1])
			corr_matrix = np.corrcoef(fit_y_pts, y_predict)
			corr = corr_matrix[0,1]
			new_rsq = corr**2
			if 1-new_rsq < 1 - rsq:
				rsq = new_rsq
				rsq_index = i
				winning_combo = comb

		# print(x_data[winning_combo])
		return fit[rsq_index]


	def calc_fractal_dimension(self,show_graph=False):
		OIsize = self.octree_size
		S0 = 1
		
		# x_IO = np.zeros((self.max_depth))
		
		fract_dim_fit = self.bestfit(self.s_data,self.Ns_data,self.bestfitlen)
		# fract_dim_fit = self.bestfit(self.s_data,self.Ns_data,self.bestfitlen,1,self.dm.radius/OIsize)

		fig, ax = plt.subplots(2,1)
		ax.flatten()


		ax[0].set_title('D = {:.2f}, unit side length'.format(fract_dim_fit[0]))
		# ax[1].set_title('D = {:.2f} for T = {}'.format(fract_dim_fit[0],self.dm.Temp))
		ax[0].plot(1/self.s_data,self.Ns_data,marker='*',label='Frac dim data')
		ax[0].loglog(1/self.s_data,np.exp(np.log(1/self.s_data)*fract_dim_fit[0]+fract_dim_fit[1]),label='log(y) = {:.2f}*log(x) + {:.2f}'.format(fract_dim_fit[0],fract_dim_fit[1]))
		ax[0].set_xlabel('log(1/(Unit side lengths))')
		ax[0].set_ylabel('log(Number of boxes to enclose)')


		ax[1].plot(1/self.s_data,np.log(self.Ns_data)/np.log(1/self.s_data))
		ax[1].set_xscale('log')
		ax[1].set_xlabel('log(1/(Unit side lengths))')
		ax[1].set_ylabel('Fractal Dimension')
		
		avg = np.mean(np.log(self.Ns_data)/np.log(1/self.s_data))
		# print('y avg: {}'.format(avg))
		ax[1].axhline(y=np.mean(np.log(self.Ns_data)/np.log(1/self.s_data)), color='g')

		ax[0].axvline(1/1, color='g')
		ax[0].axvline(1/(self.dm.radius/OIsize), color='g')

		fig.legend()
		# plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
		# mode="expand", borderaxespad=0, ncol=3)
		fig.set_figheight(10)
		fig.set_figwidth(10)
		plt.tight_layout()
		plt.savefig(self.data_folder+"FractDim.png")
		if show_graph:
			plt.show()
		plt.close()
		# fig.close()

		return fract_dim_fit[0]
	

	#function adapted from:
	#http://www.open3d.org/docs/release/python_example/geometry/octree/index.html
	def f_traverse(self,node, node_info):

		if isinstance(node, o3d.geometry.OctreePointColorLeafNode):
			self.tree_info[node_info.depth-1] += 1
		if isinstance(node, o3d.geometry.OctreeInternalPointNode):
			for child in node.children:
				if child is not None:
					self.tree_info[node_info.depth-1] += 1
					break
		#if return True, f_traverse will stop prematurely
		early_stop = False
		return early_stop


class voxelize(object):
	"""docstring for voxelize"""
	def __init__(self, data_folder,min_vox_size):
		super(voxelize, self).__init__()
		self.data_folder = data_folder
		self.dm = datamgr(self.data_folder)
		self.dm.shift_to_first_quad()
		self.min_vox_size = min_vox_size

		self.volume_elements = []
		
		self.data_range = get_data_range(self.data_folder)

		axes_len = [self.data_range[0]-self.data_range[1],self.data_range[2]-self.data_range[3],self.data_range[4]-self.data_range[5],]
		self.init_vox_size = np.max(axes_len)/2

		self.book = Tree()
		self.num_nodes = 0
		self.book.create_node(self.num_nodes,self.num_nodes)
		self.deepest_nodes = [self.book.get_node(self.num_nodes)]
		self.max_depth = self.book.depth()
		'''
		#The book implimentation stores too much redundent data
		book[chapter][sentence][word][letters]

		chapter: chapter is defined as its n value when it was written
		sentence: sentence is defined as a list of three words corresponding
					to the row, col, dep of a particular chapter.
					Basically sentence is what contains info needed to 
					determine where a box is.
					A sentence defines a colored box.
					The chapter a sentence is in defines the box size.
		word: string of 1's and 0's that define a chapter's containing boxes.
					A word can either describe a row, col, or dep.  
		letter: the 1's and 0's making up a word
		'''

		self.data_range = get_data_range(self.data_folder)

	#TODO
	def pt_to_ind(self,pt,vox_size):
		x_ind = np.ceil(pt[0]/vox_size)-1
		y_ind = np.ceil(pt[1]/vox_size)-1
		z_ind = np.ceil(pt[2]/vox_size)-1
		return (int(x_ind),int(y_ind),int(z_ind))

	def ind_to_pt(self,sentence):
		vertex = [0,0,0]
		vox_size = 0
		n=0
		for point in list(zip(*[iter(sentence)]*3)):
			for coord,letter in enumerate(point):
				vox_size = self.init_vox_size/(2**n)
				if letter == '1':
					vertex[coord] += vox_size
			n += 1

		vertices = []
		for i in range(8):
			vertices.append(vertex.copy())
			
		vertices[1][0] += vox_size
		vertices[2][1] += vox_size
		vertices[3][2] += vox_size

		vertices[4][0] += vox_size
		vertices[4][1] += vox_size
		vertices[5][1] += vox_size
		vertices[5][2] += vox_size
		vertices[6][0] += vox_size
		vertices[6][2] += vox_size
		
		vertices[7][0] += vox_size
		vertices[7][1] += vox_size
		vertices[7][2] += vox_size
		
		return vertices


	def box_contains(self,sentence,pt,radius):
		box_verts = self.ind_to_pt(sentence)
		center_to_vert_dists = []
		for vert in box_verts:
			center_to_vert_dists.append(dist(vert,pt))

		#if 1, that vertex is outside (or at) radius
		#if 0, that vertex is inside radius
		compare_radius = np.where(np.array(center_to_vert_dists)>=radius,1,0)
		#if all vertices are inside circle
		#i.e. circle completely encloses voxel
		if np.sum(compare_radius) == 0:
			# self.volume_elements.append(sentence)
			#return -1 so we can mark this node as a volume element
			return -1 #dont add to sentence since this is a volume element, dont need to see it again
		#if all vertices are outside circle
		elif np.sum(compare_radius) == compare_radius.shape:
			vox_size = self.init_vox_size/(2**self.max_depth)
			compare_vox = np.where(np.array(center_to_vert_dists)>vox_size,1,0)
			#is the circle contained in the box?
			if np.sum(compare_vox) == 0:
				return 1
			else:
				return 0
		else:
			return 1
		
	def write_chapter(self,n):
		vox_size = self.init_vox_size/(2**n)
		chapter = []

		temp_deepest_nodes = []
		for node in self.deepest_nodes:
			word = ''	
			t_node = node
			while(t_node.tag != 0):
				word = t_node.data + word
				t_node = self.book.parent(t_node.tag)

			for row in [0,1]:
				for col in [0,1]:
					for dep in [0,1]:
						new_word = str(row) + str(col) + str(dep)
						for pt in self.dm.data:
							test_ind = word + new_word
							result = self.box_contains(test_ind,pt,self.dm.radius)
							if result == 1:
								#if true, add this node
								self.num_nodes += 1
								self.book.create_node(str(self.num_nodes),str(self.num_nodes),data=new_word,parent=node)
								temp_deepest_nodes.append(self.book.get_node(str(self.num_nodes)))
								break
							elif result == -1:
								self.num_nodes += 1
								#if true, this node is a volume element
								self.book.create_node(str(self.num_nodes)+"VOL",str(self.num_nodes)+"VOL",data=new_word,parent=node)
								break


		self.deepest_nodes = temp_deepest_nodes
		self.max_depth += 1
		return vox_size

	# def 

	def write_book(self):
		curr_vox_size = self.init_vox_size
		n = 0
		while curr_vox_size > self.min_vox_size:# and n != 3:
			curr_vox_size = self.write_chapter(n)
			print(curr_vox_size,n)
			n += 1

		with open(self.data_folder+'tree.json','w') as f:
			f.write(self.book.to_json(with_data=True))
		self.book.save2file(self.data_folder+'tree.txt')
		# print(self.book)

	def rewrite_book(self,parentnode,jsonnode):
		for children in jsonnode['children']:
			for child in children:
				# print(print(children))
				# print(children[child].keys())
				self.num_nodes += 1
				new_node = self.book.create_node(str(self.num_nodes),str(self.num_nodes),parent=parentnode,data=children[child]['data'])
				if 'children' in children[child].keys():
					self.rewrite_book(new_node,children[child])

	def open_book(self):
		with open(self.data_folder+'tree.json','r') as f:
			jsonbook = json.load(f)
		
		jsonroot = jsonbook['0']
		root = self.book.get_node(self.num_nodes)
		
		self.rewrite_book(root,jsonroot)
		self.book.show()
		

	# def display_book():



# class datamgr(object):
# 	"""docstring for datamgr"""
# 	# def __init__(self, data_folder,voxel_buffer=5,ppb=3000):
# 	def __init__(self, data_folder,ppb=300000):
# 		super(datamgr, self).__init__()
# 		self.data_folder = data_folder
# 		#how many points in single ball pointcloud
# 		self.ppb = ppb
# 		self.data,self.radius,self.mass,self.moi = get_data(self.data_folder)
# 		self.nBalls = self.data.shape[0]
# 		# self.buffer = voxel_buffer # how many extra voxels in each direction 

# 	def vox_init(self,num_vox):
# 		self.data_range = get_data_range(self.data_folder)
# 		data_abs_max = max(self.data_range,key=abs) 
# 		self.vox_size = (data_abs_max*2)/num_vox
# 		# self.vox_size = 
# 		# print(self.vox_per_radius)
# 		# self.vox_rep = np.zeros((num_vox+self.buffer*2,num_vox+self.buffer*2,num_vox+self.buffer*2))

# 	#evenly spaced points on sphere code form:
# 	#http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
# 	def gen_pt_cloud(self,pt):
# 		goldenRatio = (1 + 5**0.5)/2
# 		i = np.arange(0, self.ppb)
# 		theta = 2 * np.pi * i / goldenRatio
# 		phi = np.arccos(1 - 2*(i+0.5)/self.ppb)

# 		return_array = np.zeros((self.ppb,3))
# 		return_array[:,0] = self.radius*(np.cos(theta) * np.sin(phi)) + pt[0]
# 		return_array[:,1] = self.radius*(np.sin(theta) * np.sin(phi)) + pt[1]
# 		return_array[:,2] = self.radius*np.cos(phi) + pt[2]
		
# 		return return_array

# 	def get_center_pt(self,ind):
# 		return [ind[0]*self.radius+self.radius/2,ind[1]*self.radius+self.radius/2,ind[2]*self.radius+self.radius/2]

# 	#Not tested
# 	def pt_to_ind(self,pt):
# 		x_ind = np.ceil(pt[0]/self.vox_size)-1
# 		y_ind = np.ceil(pt[1]/self.vox_size)-1
# 		z_ind = np.ceil(pt[2]/self.vox_size)-1
# 		return (int(x_ind),int(y_ind),int(z_ind))

# 	#Not tested
# 	def ind_to_pt(self,ind):
# 		vertex1 = (self.vox_size*ind[0],self.vox_size*ind[1],self.vox_size*ind[2])
# 		vertex2 = (self.vox_size*(ind[0]+1),self.vox_size*ind[1],self.vox_size*ind[2])
# 		vertex3 = (self.vox_size*ind[0],self.vox_size*(ind[1]+1),self.vox_size*ind[2])
# 		vertex4 = (self.vox_size*ind[0],self.vox_size*ind[1],self.vox_size*(ind[2]+1))
# 		vertex5 = (self.vox_size*(ind[0]+1),self.vox_size*(ind[1]+1),self.vox_size*ind[2])
# 		vertex6 = (self.vox_size*ind[0],self.vox_size*(ind[1]+1),self.vox_size*(ind[2]+1))
# 		vertex7 = (self.vox_size*(ind[0]+1),self.vox_size*ind[1],self.vox_size*(ind[2]+1))
# 		vertex8 = (self.vox_size*(ind[0]+1),self.vox_size*(ind[1]+1),self.vox_size*(ind[2]+1))

# 		return [vertex1,vertex2,vertex3,vertex4,vertex5,vertex6,vertex7,vertex8]

# 	#gotta write my own voxelization algorithm
# 	#open3d uses too much memory and has a memory leak
# 	#num_vox is number of voxels across largest dimension //TODO
# 	# def vox_me_bro(self,num_vox=500):
# 	# 	self.vox_init(num_vox)
# 	# 	self.vox_array = np.zeros((num_vox,num_vox,num_vox),dtype=np.uint8)
# 	# 	# print(self.data_range)
		
	# 	self.data[:,0] -= self.data_range[1] 
	# 	self.data[:,1] -= self.data_range[3] 
	# 	self.data[:,2] -= self.data_range[5] 


	# 	vox_in_radius = np.ceil(self.radius/self.vox_size)
	# 	print(self.radius)
	# 	print(self.vox_size)
	# 	print(vox_in_radius)
	# 	for pt in self.data:
	# 		ind = self.pt_to_ind(pt)
	# 		# box_bounds = [[ind[0]-vox_in_radius,ind[0]+vox_in_radius],[ind[1]-vox_in_radius,ind[1]+vox_in_radius],[ind[2]-vox_in_radius,ind[2]+vox_in_radius]] 
	# 		buffer = 2
	# 		box_bounds = [[ind[0]-vox_in_radius-buffer,ind[0]+vox_in_radius+buffer],[ind[1]-vox_in_radius-buffer,ind[1]+vox_in_radius+buffer],[ind[2]-vox_in_radius-buffer,ind[2]+vox_in_radius+buffer]] 
	# 		box_bounds = np.array(box_bounds,dtype=int)
	# 		print(box_bounds)
	# 		for x in range(box_bounds[0][0],box_bounds[0][1]):
	# 			for y in range(box_bounds[1][0],box_bounds[1][1]):
	# 				for z in range(box_bounds[2][0],box_bounds[2][1]):
	# 					center_pt = self.get_center_pt(ind)
	# 					print(center_pt,dist(center_pt,pt),self.radius)
	# 					for vert in self.ind_to_pt(ind):
	# 						if dist(vert,pt) < self.radius and self.vox_array[x,y,z] == 0:
	# 							self.vox_array[x,y,z] = 1
	# 							break

	# 		print(np.sum(self.vox_array))
	# 		print(ind)
	# 		exit(0)
	# 		self.vox_array[self.pt_to_ind(pt)] += 1

		# np.fromiter((self))


		# print("HERE")
		# ax = plt.figure().add_subplot(projection='3d')
		# print("HERE1")
		# ax.voxels(self.vox_array, edgecolor='k')
		# print("HERE2")

		# plt.show()



	#adapted from
	#https://towardsdatascience.com/how-to-voxelize-meshes-and-point-clouds-in-python-ca94d403f81d
	# def vox_me_bro(self,num_vox=5000):
	# 	num_vox = 5000

	# 	# self.point_cloud = np.zeros((2*self.ppb,3))
	# 	self.point_cloud = np.zeros((self.nBalls*self.ppb,3))
	# 	# self.point_cloud = np.zeros((self.ppb,3))
	# 	self.vox_init(num_vox)

	# 	# colors = np.random.rand((self.ppb,3))	
	# 	pcd = o3d.geometry.PointCloud()
	# 	for ind,pt in enumerate(self.data):
	# 		# self.gen_pt_cloud(pt)
	# 		self.point_cloud[ind*self.ppb:(ind+1)*self.ppb] = self.gen_pt_cloud(pt)
	# 		# Add the points, colors and normals as Vectors
	# 	# pcd.points = o3d.utility.Vector3dVector(self.gen_pt_cloud(pt))
	# 	pcd.points = o3d.utility.Vector3dVector(self.point_cloud)
	# 	print(pcd)
	# 	# exit(0)
	# 	pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=(self.ppb*self.nBalls, 3)))
		
	# 	print("HERE1")

	# 	# o3d.visualization.draw_geometries([pcd])
	# 	# tetra_mesh, pt_map = o3d.geometry.TetraMesh.create_from_point_cloud(pcd)
	# 	print("HERE2")
	# 	alpha = 0.03
	# 	# mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
	# 	# mesh.compute_vertex_normals()
	# 	# o3d.visualization.draw_geometries([mesh], mesh_show_back_face=True)
	# 	# mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
	# 	# mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha,tetra_mesh,pt_map)
	# 	print("HERE3")
	# 	# mesh.compute_vertex_normals()

	# 	# self.point_cloud[ind*self.ppb:(ind+1)*self.ppb] = self.gen_pt_cloud(pt)
	# 	voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd,voxel_size=self.vox_size)
	# 	print(voxel_grid)
	# 	# voxel_grid = o3d.geometry.VoxelGrid.create_dense(voxel_grid)
	# 	# print(voxel_grid.get_voxels)
		

	# 	vis = o3d.visualization.Visualizer()
	# 	# # Create a window, name it and scale it
	# 	vis.create_window(window_name='VoxVis', width=800, height=600)
	# 	# vis.draw_geometries([mesh], mesh_show_back_face=True)

	# 	# # # Add the voxel grid to the visualizer
	# 	vis.add_geometry(voxel_grid)

	# 	# # # We run the visualizater
	# 	vis.run()
	# 	# # # Once the visualizer is closed destroy the window and clean up
	# 	vis.destroy_window()
	# 	exit(0)


	#adapted from
	#https://towardsdatascience.com/how-to-voxelize-meshes-and-point-clouds-in-python-ca94d403f81d
	# def vox_me_bro1(self,vox_res=50,cubic_size=1.0):
	# 	self.point_cloud = np.zeros((self.nBalls*self.ppb,3))
	# 	# self.point_cloud = np.zeros((self.ppb,3))
	# 	# self.vox_init(num_vox)

	# 	# colors = np.random.rand((self.ppb,3))	

	# 	#check if pointcloud exists
	# 	vox_file = self.data_folder + "voxFile.out"
	# 	# vox_file = self.data_folder + "voxFile.pts"
	# 	# if os.path.isfile(vox_file):
	# 	if False:
	# 		# pcd = o3d.io.read_point_cloud(vox_file)
	# 		loaded_arr = np.loadtxt(vox_file) 
	# 		voxel_array = loaded_arr.reshape(\
 #    			loaded_arr.shape[0], loaded_arr.shape[1] // vox_res,vox_res)
	# 	else:
	# 		# for ind,pt in enumerate([self.data[0]]):
	# 		for ind,pt in enumerate(self.data):
	# 			# self.gen_pt_cloud(pt)
	# 			self.point_cloud[ind*self.ppb:(1+ind)*self.ppb] = self.gen_pt_cloud(pt)
	# 			# Add the points, colors and normals as Vectors
	# 		pcd = o3d.geometry.PointCloud()
	# 		pcd.points = o3d.utility.Vector3dVector(self.point_cloud)
	# 		# pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=(self.ppb*self.nBalls, 3)))
	# 		visualization = True
	# 		cubic_size = cubic_size
	# 		# voxel_resolution = 50.0
	# 		voxel_resolution = vox_res

	# 		voxel_grid, voxel_carving, voxel_surface = vox_carve(pcd, cubic_size, voxel_resolution)

	# 		# # pcd.points = o3d.utility.Vector3dVector(self.gen_pt_cloud(pt))
	# 		# pcd.estimate_normals()
	# 		# # exit(0)
		
	# 		# print("HERE1")
	# 		# o3d.io.write_voxel_grid(vox_file, voxel_grid)
	# 		# print(voxel_grid.get_geometry_type())
	# 		print(len(voxel_grid.get_voxels()))

	# 		voxel_array = np.zeros((int(voxel_resolution),int(voxel_resolution),int(voxel_resolution)))
	# 		voxels = voxel_grid.get_voxels()

	# 		for vox in voxels:
	# 			vxpt = vox.grid_index
	# 			voxel_array[vxpt[0],vxpt[1],vxpt[2]] += 1

	# 		print(voxel_array)
	# 		print(np.sum(voxel_array))

	# 		# np.savetxt(vox_file,voxel_array.reshape(voxel_array.shape[0],-1))

	# 	# o3d.visualization.draw_geometries([pcd])
	# 	# tetra_mesh, pt_map = o3d.geometry.TetraMesh.create_from_point_cloud(pcd)
	# 	# print("HERE2")
	# 	# alpha = 0.03
	# 	# mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
	# 	# mesh.compute_vertex_normals()
	# 	# o3d.visualization.draw_geometries([mesh], mesh_show_back_face=True)
	# 	# # mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
	# 	# # mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha,tetra_mesh,pt_map)
	# 	# print("HERE3")
	# 	# exit(0)
	# 	# mesh.compute_vertex_normals()

		

	# 	# self.point_cloud[ind*self.ppb:(ind+1)*self.ppb] = self.gen_pt_cloud(pt)
	# 	# voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd,voxel_size=0.000001)
	# 	# voxel_grid = o3d.geometry.VoxelGrid.create_dense(voxel_grid)
	# 	# print(voxel_grid.get_voxels)
		





	# 	o3d.visualization.draw_geometries([voxel_grid])
	# 	# vis = o3d.visualization.Visualizer()
	# 	# # # Create a window, name it and scale it
	# 	# vis.create_window(window_name='VoxVis', width=800, height=600)
	# 	# # vis.draw_geometries([mesh], mesh_show_back_face=True)
	# 	# o3d.visualization.draw_geometries([voxel_carving])
	# 	# # # # Add the voxel grid to the visualizer
	# 	# vis.add_geometry(voxel_grid)

	# 	# # # # We run the visualizater
	# 	# vis.run()
	# 	# # # # Once the visualizer is closed destroy the window and clean up
	# 	# vis.destroy_window()



	# 	# ax = plt.figure().add_subplot(projection='3d')
	# 	# ax.voxels(voxel_array, edgecolor='k')

	# 	plt.show()
		

	# 	exit(0)

				
		
	
