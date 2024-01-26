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
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u




folder1 = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_data/jobs/errorckcsvlognorm1/N_10/T_3/"
folder2 = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_data/jobs/errorckh5lognorm1/N_10/T_3/"

data1,radius1,mass1,moi1 = u.get_data(folder1,9)
data2,radius2,mass2,moi2 = u.get_data(folder2,9)

print(mass1)
print(mass2)

