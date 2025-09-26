import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import exists

def read_single_line(filename, line_number):
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i == line_number:
                values = line.strip().split(',')
                float_array = np.array([float(v) for v in values])
                return float_array
    raise ValueError(f"Line {line_number} not found in {filename}")

def calc_COM(pos,mass):
	Mtot = np.sum(mass,dtype=np.float128)
	if len(pos.shape) == 3:
		COM = []
		for i in range(pos.shape[0]):
			mx = np.sum(np.multiply(pos[i,:,0],mass))
			my = np.sum(np.multiply(pos[i,:,1],mass))
			mz = np.sum(np.multiply(pos[i,:,2],mass))
			COM.append(np.array([mx,my,mz],dtype=np.float128)/Mtot)
		return np.array(COM)
	elif len(pos.shape) == 2:
		mx = np.sum(np.multiply(pos[:,0],mass))
		my = np.sum(np.multiply(pos[:,1],mass))
		mz = np.sum(np.multiply(pos[:,2],mass))
		return np.array(np.array([mx,my,mz],dtype=np.float128)/Mtot)
	else:
		raise ValueError("ERROR: pos should be 2 or 3 dimensional but has a shape of {}".format(pos.shape))

def tot_omega(pos,vel,mass):
	for i in range(pos.shape[0]):
		COM = calc_COM(pos[i],mass)
		for j in range(pos.shape[1]):
			print('Ball {}: {}'.format(j,mass[j]*np.cross((COM-pos[i,j,:]),vel[i,j,:])))
		if i > 10:
			exit(0)
def main():
	# new_data = True
	get_torque = False
	get_slide = False
	get_roll = False
	get_pos = True
	get_vel = True
	get_mass = True

	path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/step_back10/N_2/T_1/"
	sims = [1]
	for j,sim in enumerate(sims):
		if sim == 0:
			fileprefix = "2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
		else:
			fileprefix = "{}_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_".format(sim)

		# sav_torque = path + 'torque.csv'
		# sav_slide = path + 'sliding.csv'
		# sav_roll = path + 'rolling.csv'

		# has_torque = exists(sav_torque)
		# has_slid = exists(sav_slide)
		# has_roll = exists(sav_roll)
		if True:
			if get_mass:
				massfile = path+fileprefix + "constants.csv"
				massdata = np.loadtxt(massfile,delimiter=',',dtype=np.float128)
				mass = massdata[:,1]

			if get_pos:
				posfile = path + fileprefix + "posData.csv"
				posdata = np.loadtxt(posfile,delimiter=',',dtype=np.float128,skiprows=1)
				pos = posdata.reshape(posdata.shape[0],int(posdata.shape[1]/3),3)

			if get_vel:
				velfile = path + fileprefix + "velData.csv"
				veldata = np.loadtxt(velfile,delimiter=',',dtype=np.float128,skiprows=1)
				vel = veldata.reshape(posdata.shape[0],int(posdata.shape[1]/3),3)

			if get_torque:
				# if new_data or (not has_slid or not has_roll or not has_torque):
				torquefile = path + fileprefix + "torqueData.csv"
				torquedata = np.loadtxt(torquefile,delimiter=',',dtype=np.float64,skiprows=1)
				torquedata = torquedata.reshape(torquedata.shape[0],int(torquedata.shape[1]/3),3)
				torque = np.zeros(pos.shape)
				
			if get_slide:
				slidefile = path + fileprefix + "slideData.csv"
				slidedata = np.loadtxt(slidefile,delimiter=',',dtype=np.float64,skiprows=1)
				slidedata = slidedata.reshape(slidedata.shape[0],int(slidedata.shape[1]/3),3)
				slide = np.zeros(pos.shape)
				
			if get_roll:
				rollfile = path + fileprefix + "rollData.csv"
				rolldata = np.loadtxt(rollfile,delimiter=',',dtype=np.float64,skiprows=1)
				rolldata = rolldata.reshape(rolldata.shape[0],int(rolldata.shape[1]/3),3)
				roll = np.zeros(pos.shape)


			for i in range(0,pos.shape[1]):
				if get_torque:
					torque[:,i,:] = np.sum(torquedata[:,i*pos.shape[1]:(i+1)*pos.shape[1],:],1)
				if get_slide:
					slide[:,i,:] = np.sum(slidedata[:,i*pos.shape[1]:(i+1)*pos.shape[1],:],1)
				if get_roll:	
					roll[:,i,:] = np.sum(rolldata[:,i*pos.shape[1]:(i+1)*pos.shape[1],:],1)

			# np.savetxt(sav_torque,torque.reshape(),delimiter=',')
			# np.savetxt(sav_slide,slide,delimiter=',')
			# np.savetxt(sav_roll,roll,delimiter=',')
		else:
			torque = np.loadtxt(sav_torque,delimiter=',')
			slide = np.loadtxt(sav_slide,delimiter=',')
			roll = np.loadtxt(sav_roll,delimiter=',')


		# tot_omega(pos,vel,mass)
		print(calc_COM(pos[100000:101000],mass))

		# i=1
		# j=0
		# print(torque[i,j,:])
		# print(np.cross(pos[i,j,:],(slide[i,j,:]+roll[i,j,:])))
		# print(torque[i,j,:] - np.cross(pos[i,j,:],(slide[i,j,:]+roll[i,j,:])))
		# print(np.sum(torque[i,:,:] - np.cross(pos[i,:,:],(slide[i,:,:]+roll[i,:,:])),0))


if __name__ == '__main__':
	main()