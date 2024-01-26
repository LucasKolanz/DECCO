import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
import csv

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
sys.path.append(project_path)
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u
import porosity_FD as p
# import utils_old as u



def main():
	folder = "/home/lucas/Desktop/SpaceLab_data/jobs/overlapError2/N_5/T_1000/"

	print(p.number_of_contacts(folder,4))

if __name__ == '__main__':
	main()




