import numpy as np
import argparse
import os
import time
import sys
file_path = os.path.abspath(__file__)  # Get absolute path of the file
folder = os.path.dirname(file_path)  # Get the directory.

project_path = '/'.join(folder.split('/')[0:folder.split('/').index("SpaceLab_data")])+'/SpaceLab/'


sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

density = 2.25 #g/cm^3

def dataOut(x,y,z,r,folder=""):
	if folder == "":
		file_path = os.path.abspath(__file__)  # Get absolute path of the file
		folder = os.path.dirname(file_path)+'/'  # Get the directory
	numParticles = len(x)
	mass = density*(4/3)*np.pi*r**3
	moi = (2/5)*mass*r*r
	u.write_constants(folder+f"{numParticles}_constants.csv",r,mass,moi)

	u.write_energy(folder+f"{numParticles}_energy.csv",0,0,0,0,0,0)

	pos = np.zeros(shape=(numParticles,3))
	w = np.copy(pos)
	vel = w
	pos[:,0] = x
	pos[:,1] = y
	pos[:,2] = z

	u.write_simData(folder+f"{numParticles}_simData.csv",pos,w,vel)

def rCM(x,y,z,r):
	xCM=np.sum(x*r**3)/np.sum(r**3)
	yCM=np.sum(y*r**3)/np.sum(r**3)
	zCM=np.sum(z*r**3)/np.sum(r**3)
	return(xCM,yCM,zCM) # this just calculates the center of mass


def constant(numParticles, radius):
	# the first particle is set at the origin and has radius 1
	# for now, all radii are 1
	x=np.array([0])
	y=np.array([0])
	z=np.array([0])
	r=np.array([radius])
	nTot = numParticles-1
	# now we add the particles one by one
	for i in range(nTot):
	# first we find the radius of a sphere that completely contains the current cluster.
	# remember the cluster is set in such a way that the center of mass is always in the origin
	# note the factor 2.1 because the guess particle radius needs also to be considered
		dist=np.sqrt((x)**2+(y)**2+(z)**2)
		rArea=dist.max()+2.1*r.max()
	# now we simulate a large number of particles with random velocities on the surface of that sphere
	# first we simulate a filled cube
		nSim=100
		xSim=(np.random.rand(nSim)-0.5)*rArea*2
		ySim=(np.random.rand(nSim)-0.5)*rArea*2
		zSim=(np.random.rand(nSim)-0.5)*rArea*2
	#The we select those inside a sphere
		jj=np.where((xSim**2+ySim**2+zSim**2)<=(rArea**2))
		xSim=xSim[jj]; ySim=ySim[jj]; zSim=zSim[jj]
	# finally we put them all on the sphere
		rr=np.sqrt(xSim**2+ySim**2+zSim**2)
		xSim=xSim/rr*rArea ; ySim=ySim/rr*rArea ; zSim=zSim/rr*rArea
	# now we do the same for their velocities
		vxSim=(np.random.rand(nSim)-0.5)*2
		vySim=(np.random.rand(nSim)-0.5)*2
		vzSim=(np.random.rand(nSim)-0.5)*2
		jj=np.where((vxSim**2+vySim**2+vzSim**2)<=(1))
		vxSim=vxSim[jj]; vySim=vySim[jj]; vzSim=vzSim[jj]
		vv=np.sqrt(vxSim**2+vySim**2+vzSim**2)
		vxSim=vxSim/vv ; vySim=vySim/vv ; vzSim=vzSim/vv
	# now we make sure we have as many velocities as we have positions
		jj=np.minimum(xSim.size,vxSim.size)
		xSim=xSim[:jj] ; ySim=ySim[:jj] ; zSim=zSim[:jj]
		vxSim=vxSim[:jj] ; vySim=vySim[:jj] ; vzSim=vzSim[:jj]
	# here we give all the guess particles the radius 1
		rSim=xSim*0+radius
		
	# here we check which of the guess particles hits. Once we find one we stop 
	# checking and we see where that particle will position itself.
		for ig in range(xSim.size):
			posSim=np.array([xSim[ig],ySim[ig],zSim[ig]])
			velSim=np.array([vxSim[ig],vySim[ig],vzSim[ig]])
			dists=np.zeros(x.size)
			for ip in range(x.size):
				posPart=np.array([x[ip],y[ip],z[ip]])
				dists[ip]=np.linalg.norm(np.cross(velSim,posSim-posPart))
			kk=np.where(dists<(r+rSim[ig]))
			if (kk[0].size>0): break
	# first we check if the velocity should be reversed.
		if np.dot(posSim,velSim)>0: velSim=-velSim
	# we now create a trajectory to see where the particle will land.
		# dt=0.001
		dt=radius*0.1
		dists=np.sqrt((posSim[0]-x)**2+(posSim[1]-y)**2+(posSim[2]-z)**2)
		jj=np.where(dists<(r+rSim[ig]))
		while jj[0].size==0:
			posSim+=velSim*dt
			dists=np.sqrt((posSim[0]-x)**2+(posSim[1]-y)**2+(posSim[2]-z)**2)
			jj=np.where(dists<(r+rSim[ig]))
	# we have now found the position of the new particle where it sticks
		x=np.concatenate([x,[posSim[0]]])
		y=np.concatenate([y,[posSim[1]]])
		z=np.concatenate([z,[posSim[2]]])
		r=np.concatenate([r,[rSim[ig]]])
		xCM,yCM,zCM=rCM(x,y,z,r)
		x=x-xCM
		y=y-yCM
		z=z-zCM
		print(x.size)

	dataOut(x,y,z,r)

def lognormal(numParticles, radius, sigma):

	

	nTot=numParticles-1 # the total number of particles in the final product

	#We want the average radius to still be 1e-5 cm (or whatever -r is set to)
	a_max = radius*np.exp(-5*sigma**2/2)
	#properties of the Gaussian used for the distribution
	meanG=np.log(a_max)-sigma**2 #
	sigG=sigma  # 0.2 is a reasonable value

	# the first particle is set at the origin and has radius 1
	# for now, all radii are 1
	x=np.array([0])
	y=np.array([0])
	z=np.array([0])
	r=np.array(np.random.lognormal(size=1,mean=meanG,sigma=sigG))

	# now we add the particles one by one
	for i in range(nTot):
		# first we find the radius of a sphere that completely contains the current cluster.
		# remember the cluster is set in such a way that the center of mass is always in the origin
		# note the factor 2.1 because the guess particle radius needs also to be considered
		dist=np.sqrt((x)**2+(y)**2+(z)**2)
		rArea=dist.max()+2.1*r.max()
		# now we simulate a large number of particles with random velocities on the surface of that sphere
		# first we simulate a filled cube
		nSim=100
		xSim=(np.random.rand(nSim)-0.5)*rArea*2
		ySim=(np.random.rand(nSim)-0.5)*rArea*2
		zSim=(np.random.rand(nSim)-0.5)*rArea*2
		#The we select those inside a sphere
		jj=np.where((xSim**2+ySim**2+zSim**2)<=(rArea**2))
		xSim=xSim[jj]; ySim=ySim[jj]; zSim=zSim[jj]
		# finally we put them all on the sphere
		rr=np.sqrt(xSim**2+ySim**2+zSim**2)
		xSim=xSim/rr*rArea ; ySim=ySim/rr*rArea ; zSim=zSim/rr*rArea
		# now we do the same for their velocities
		vxSim=(np.random.rand(nSim)-0.5)*2
		vySim=(np.random.rand(nSim)-0.5)*2
		vzSim=(np.random.rand(nSim)-0.5)*2
		jj=np.where((vxSim**2+vySim**2+vzSim**2)<=(1))
		vxSim=vxSim[jj]; vySim=vySim[jj]; vzSim=vzSim[jj]
		vv=np.sqrt(vxSim**2+vySim**2+vzSim**2)
		vxSim=vxSim/vv ; vySim=vySim/vv ; vzSim=vzSim/vv
		# now we make sure we have as many velocities as we have positions
		jj=np.minimum(xSim.size,vxSim.size)
		xSim=xSim[:jj] ; ySim=ySim[:jj] ; zSim=zSim[:jj]
		vxSim=vxSim[:jj] ; vySim=vySim[:jj] ; vzSim=vzSim[:jj]

		rSim=np.random.lognormal(size=xSim.size,mean=meanG,sigma=sigG)
		
		# here we check which of the guess particles hits. Once we find one we stop 
		# checking and we see where that particle will position itself.
		for ig in range(xSim.size):
			posSim=np.array([xSim[ig],ySim[ig],zSim[ig]])
			velSim=np.array([vxSim[ig],vySim[ig],vzSim[ig]])
			dists=np.zeros(x.size)
			for ip in range(x.size):
				posPart=np.array([x[ip],y[ip],z[ip]])
				dists[ip]=np.linalg.norm(np.cross(velSim,posSim-posPart))
			kk=np.where(dists<(r+rSim[ig]))
			if (kk[0].size>0): 
				break
		# first we check if the velocity should be reversed.
		if np.dot(posSim,velSim)>0: 
			velSim=-velSim
		# we now create a trajectory to see where the particle will land.
		# dt=0.001
		dt=0.1*np.min(np.array([np.min(r),rSim[ig]]))
		dists=np.sqrt((posSim[0]-x)**2+(posSim[1]-y)**2+(posSim[2]-z)**2)
		jj=np.where(dists<(r+rSim[ig]))
		while jj[0].size==0:
			posSim+=velSim*dt
			dists=np.sqrt((posSim[0]-x)**2+(posSim[1]-y)**2+(posSim[2]-z)**2)
			jj=np.where(dists<(r+rSim[ig]))
		# we have now found the position of the new particle where it sticks
		x=np.concatenate([x,[posSim[0]]])
		y=np.concatenate([y,[posSim[1]]])
		z=np.concatenate([z,[posSim[2]]])
		r=np.concatenate([r,[rSim[ig]]])
		xCM,yCM,zCM=rCM(x,y,z,r)
		x=x-xCM
		y=y-yCM
		z=z-zCM
		print(x.size)

	dataOut(x,y,z,r)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--numParticles', default=300, help="The number of particles in the final aggregate.")
	parser.add_argument('-r', '--radius', default=1e-5, help="The radius (in cm) of particles (only for constant distribution).")
	parser.add_argument('-s', '--sigma', default=0.2, help="The sigma for the lognormal distribution (only for lognormal distribution).")
	# parser.add_argument('-m', '--mean', default=0.0, help="The mean of the lognormal distribution (this is actually the mean of the gaussian in the lognormal distribution) (only for lognormal distribution).")
	parser.add_argument('-d', '--distribution', default="lognormal", help="The mean of the lognormal distribution (this is actually the mean of the gaussian in the lognormal distribution) (only for lognormal distribution).")
	parser.add_argument('-S', '--seed', default="", help="Seed for numpy.random. Defaults to the current time.")


	args = parser.parse_args()

	if args.seed == "":
		np.random.seed(int(time.time()))
	else:
		np.random.seed(int(args.seed))


	if args.distribution == "constant":
		constant(int(args.numParticles), float(args.radius))
	elif args.distribution == "lognormal":
		lognormal(int(args.numParticles), float(args.radius), float(args.sigma))
	else:
		print("ERROR: distribution not recognized")

	file_path = os.path.abspath(__file__)  # Get absolute path of the file
	folder = os.path.dirname(file_path)  # Get the directory
	os.system(f"touch {folder}/{args.numParticles}_checkpoint.txt")
	os.system(f"touch {folder}/timing.txt")