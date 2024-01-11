import numpy as np


def main():
	start_agg_path = "/global/homes/l/lpkolanz/SpaceLab/jobs/collidable_aggregate/"
	start_agg_filebase = "_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	start_agg_size = 299

	copy_to_path = "/global/u2/l/lpkolanz/SpaceLab/jobs/weakScaleGrow1/thread_16/"
	desired_agg_size = 300
	parameters = 11

	simData = np.loadtxt(start_agg_path+str(start_agg_size)+start_agg_filebase+'simData.csv',skiprows=1,delimiter=',',dtype=np.float64)[-1]
	newSimData = simData[:parameters*(desired_agg_size+2)]
	np.savetxt(copy_to_path+str(desired_agg_size-1)+start_agg_filebase+"simData.csv",newSimData,newline=',')

	constants = np.loadtxt(start_agg_path+str(start_agg_size)+start_agg_filebase+'constants.csv',skiprows=1,delimiter=',')
	# print(constants)
	newConstants = constants[:desired_agg_size+2]
	np.savetxt(copy_to_path+str(desired_agg_size-1)+start_agg_filebase+"constants.csv",newConstants,delimiter=',')
	

if __name__ == '__main__':
	main()