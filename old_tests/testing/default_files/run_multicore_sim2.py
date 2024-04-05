import os
import sys
import subprocess


def main():
	location = sys.argv[1]
	if len(sys.argv) > 2:
		num_balls = sys.argv[2]
	else:
		num_balls = 100
	# os.system("cd {}".format(location))
	# curr_folder = os.getcwd() + '/'
	output_file = location + "sim_output.txt"
	error_file = location + "sim_errors.txt"
	# out = os.system("./ColliderSingleCore.o {}".format(curr_folder))
	# out = os.system("./ColliderSingleCore.o {} 1>> {} 2>> {}".format(curr_folder,output_file,error_file))
	
	# cmd = ["srun","-n","1","-c","2","{}ColliderSingleCore.x".format(location), location, str(num_balls)]
	cmd = ["{}ColliderMultiCore2.x".format(location), location, str(num_balls)]
	with open(output_file,"a") as out, open(error_file,"a") as err:
		subprocess.run(cmd,stdout=out,stderr=err)


if __name__ == '__main__':
	main()