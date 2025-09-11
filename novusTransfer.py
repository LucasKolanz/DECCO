"""
This file was originally written for SpaceLab/DECCO to trasnfer jobs from the novus cluster to my local machine

Author: Lucas Kolanz

This file transfers all (finished) simulations that match a specified pattern from the Novus cluster to the local computer
this file is run from. The pattern of the folders will be kept the same, but will go into a parent directory in data_directory called "jobsNovus/" 
It does this with scp and paramiko python libraries. Note that this file assumes you have an ssh config file that
specifies the hostname, user, and the location of the identity file you use to ssh into the (in this case COSINE) system.

This file is basically the same as cosineTransfer.py except that the folder patterns are slightly different.

TODO: combine this and cosineTransfer into one file that takes a command line arg for which cluster to grab from

"""

from paramiko import SSHClient, SSHConfig, AutoAddPolicy
from scp import SCPClient, SCPException
import os
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

#IF files is empty then just copy remote_file_path, otherwise loop through and copy all files from remote_file_path
def scp_transfer_with_ssh_config(hostname, remote_file_path, local_path,recursive=False): 
	# Load SSH config
	ssh_config = SSHConfig()
	user_config_file = os.path.expanduser("~/.ssh/config")
	if os.path.exists(user_config_file):
		with open(user_config_file) as f:
			ssh_config.parse(f)

	user_config = ssh_config.lookup(hostname)

	# Establish SSH connection
	ssh = SSHClient()
	ssh.set_missing_host_key_policy(AutoAddPolicy())
	ssh.connect(
		hostname=user_config.get('hostname', hostname),
		username=user_config.get('user'),
		key_filename=user_config.get('identityfile'),
		port=user_config.get('port', 22)
	)

	file_found = False

	# SCP transfer
	with SCPClient(ssh.get_transport()) as scp:
		try:
			scp.get(remote_file_path, local_path,recursive=recursive)

			file_found = True
		except SCPException as e:
			# Check if the exception is due to a missing file
			if "No such file or directory" in str(e):
				if remote_file_path[-10:] == "timing.txt":
					print(f"Sim not finished: {remote_file_path}")
				else:
					print(f"file doesn't exit: {remote_file_path}")
			else:
				# Re-raise the exception if it's not a file not found error
				raise


	# Close SSH connection
	ssh.close()

	return file_found



def list_remote_files(hostname, remote_directory):
	# Load SSH config
	ssh_config = SSHConfig()
	user_config_file = os.path.expanduser("~/.ssh/config")
	if os.path.exists(user_config_file):
		with open(user_config_file) as f:
			ssh_config.parse(f)

	user_config = ssh_config.lookup(hostname)

	# Establish SSH connection
	ssh = SSHClient()
	ssh.set_missing_host_key_policy(AutoAddPolicy())
	ssh.connect(
		hostname=user_config.get('hostname', hostname),
		username=user_config.get('user'),
		key_filename=user_config.get('identityfile'),
		port=user_config.get('port', 22)
	)

	# Command to list files in the remote directory
	command = f"ls {remote_directory}"
	stdin, stdout, stderr = ssh.exec_command(command)

	# Read the output (list of files)
	file_list = stdout.readlines()

	# Close the SSH connection
	ssh.close()

	# Process and return the file list
	return [filename.strip('\n') for filename in file_list]



def main():


	# curr_folder = os.getcwd() + '/'

	
	remote_base_folder = '/home/kolanzl/novus/kolanzl/SpaceLab_data/'

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	data_directory = input_json["data_directory"]

	job_set_names = ["const","lognorm"]
	job_set_names = ["const"]
	job_set_names = ["lognorm"]
	job_set_names = ["BAPA"]
	

	attempts = [i for i in range(30)]
	# attempts = [0,1]
	attempts_300 = attempts

	jobfolder = ""


	N = [300]
	Temps = [1000]
	M = [3,5,10,15]

	# Temps = [3]
	for j_i,job_set_name in enumerate(job_set_names):
		if job_set_name == "lognorm":
			jobfolder = "jobsCosine"
		elif job_set_name == "const":
			jobfolder = "jobsNovus"
		elif job_set_name == "BAPA":
			jobfolder = "jobs"
		else:
			print("ERROR: unrecognized job_set_name.")
			exit(0)
		for n in N:
			for m in M:
				for Temp in Temps:
					for attempt in attempts:
						
						local_job_folder = data_directory  + jobfolder + '/' + job_set_name + '_' + str(attempt)\
									+ "/M_" + str(m) + '/N_' + str(n) + '/T_' + str(Temp) + '/'
						remote_job_folder = remote_base_folder + 'jobs/' + job_set_name + '_' + str(attempt)\
									+ "/M_" + str(m) + '/N_' + str(n) + '/T_' + str(Temp) + '/'


						if os.path.exists(local_job_folder+"timing.txt"): #Have we already copied this job over?
							print(f"Job already copied: {remote_job_folder}")
							continue

						#If you want to scp a whole folder recursivly it will copy the final folder 
						#which is stupid. To counter this, take out the last folder in local_path
						#since these are the same folder
						local_job_folder = '/'.join(local_job_folder.split('/')[:-2]) + '/'
						

						# remote_files_exists = list_remote_files('Novus', remote_job_folder)
						# if len(remote_files_exists) > 1:
						# 	print(f"Remote folder {remote_job_folder} has {len(remote_files_exists)} files")

						remote_file_exists = list_remote_files('Novus', remote_job_folder+"timing.txt")
						if len(remote_file_exists) > 0 and "timing.txt" == remote_file_exists[0].split('/')[-1]:
							if not os.path.exists(local_job_folder): #folder doesnt exist locally so make the folder(s)
								os.makedirs(local_job_folder)
							print(f"Copying {remote_job_folder}")
							scp_transfer_with_ssh_config('Novus',remote_job_folder,local_job_folder,recursive=True)



if __name__ == '__main__':
	main()