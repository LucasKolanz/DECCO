import numpy as np
# import matplotlib.pyplot as plt
import os 
import sys
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u


def check_halves(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        total_lines = len(lines)

        # Check if the total number of lines is odd
        if total_lines % 2 != 0:
            return False

        half_lines = total_lines // 2

        first_half = lines[:half_lines]
        second_half = lines[-half_lines:]

        return first_half == second_half

    except FileNotFoundError:
        print(f"The file {filename} does not exist.")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False

def replace_file(filename, N):
    # Ensure N is positive
    if N <= 0:
        raise ValueError("N must be a positive integer")

    try:
        # First Pass: Count the total number of lines in the file
        with open(filename, 'r') as file:
            total_lines = sum(1 for _ in file)

        # Only proceed if the file has more than N lines
        if total_lines > N and check_halves(filename):

            # Second Pass: Keep only the last N lines
            otherfile = filename.replace("RELAX","")
            command = f"cp {otherfile} {filename}"
            # command = filename.replace("RELAX","")
            print(f"{total_lines}: {command}")
            # exit(0)
            # os.system(command)

    # except FileNotFoundError:
    #     print(f"The file {filename} does not exist.")
    except Exception as e:
    	pass
    #     print(f"An error occurred: {e}")


def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsCosine/lognorm_relax'

	dataset_name = data_prefolder.split("/")[-1]

	# sav = path+'data/{}_averageData.csv'.format(dataset_name)
	# # figure_folder = 'figuresCompare/'
	# figure_folder = path+'data/figures/'


	temps = [3,10,30,100,300,1000]
	# temps = [1000]
	Nums = [30,100,300]
	Nums = [300]
	
	
	attempts = [i for i in range(30)]

	dir_pattern = data_prefolder+"$a$/N_$n$/T_$t$/"

	for a,attempt in enumerate(attempts):
		for n,N in enumerate(Nums):
			for t,temp in enumerate(temps):
				data_folder = dir_pattern.replace("$a$",str(attempt)).replace("$n$",str(N)).replace("$t$",str(temp))

				replace_file(data_folder+"297_RELAXconstants.csv",300)



if __name__ == '__main__':
	main()