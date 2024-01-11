import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kolanzl/Desktop/SpaceLab")
import utils as u

# def make_vox_rep(data_folder,vox_size):


# def voxelize(data_folder,min_vox_size):
# 	dmgr = u.datamgr(data_folder)
# 	dmgr.vox_me_bro(num_vox)

# 	exit(0)
	

if __name__ == '__main__':
	o3do = u.o3doctree()
	o3do.test_menger_sponge()