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

	
	# novus_base_folder = '/home/kolanzl/novus/kolanzl/SpaceLab_data/'
	cosine_data_directory = '/home/physics/kolanzl/SpaceLab_data/'

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	novus_data_directory = input_json["data_directory"]

	job_set_name = "lognorm"

	attempts = [i for i in range(30)]
	attempts = [0]
	attempts_300 = attempts


	N = [30,100,300]
	N=[300]

	Temps = [3,10,30,100,300,1000]
	Temps = [3]

	for n in N:
		for Temp in Temps:
			temp_attempt = attempts
			if n == 300:
				temp_attempt = attempts_300
			for attempt in temp_attempt:
				novus_job_directory = novus_data_directory  + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				cosine_job_directory = cosine_data_directory + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

				if os.path.exists(novus_job_directory+"timing.txt"): #Have we already copied this job over?
					print(f"Job already copied: {remote_job_folder}")
					continue

				#If you want to scp a whole folder recursivly it will copy the final folder 
				#which is stupid. To counter this, take out the last folder in local_path
				#since these are the same folder
				novus_job_directory = '/'.join(novus_job_directory.split('/')[:-2]) + '/'
				

				remote_file_exists = list_remote_files('Cosine', cosine_job_directory+"timing.txt")
				if len(remote_file_exists) > 0 and "timing.txt" == remote_file_exists[0].split('/')[-1]:
					if not os.path.exists(novus_job_directory): #folder doesnt exist locally so make the folder(s)
						os.makedirs(novus_job_directory)
					print(f"Copying {cosine_job_directory} to {novus_job_directory}")
					# scp_transfer_with_ssh_config('Novus',cosine_job_directory,novus_job_directory,recursive=True)


if __name__ == '__main__':
	main()