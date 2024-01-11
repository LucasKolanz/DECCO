from paramiko import SSHClient, SSHConfig, AutoAddPolicy
from scp import SCPClient, SCPException
import os

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


	curr_folder = os.getcwd() + '/'
	remote_base_folder = '/home/physics/kolanzl/SpaceLab/'

	job_set_name = "lognorm"

	attempts = [i for i in range(30)]
	# attempts = [0,1]
	attempts_300 = attempts


	N = [30,100,300]
	# N=[100]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	error_folders = []
	for n in N:
		for Temp in Temps:
			temp_attempt = attempts
			if n == 300:
				temp_attempt = attempts_300
			for attempt in temp_attempt:
				local_job_folder = curr_folder + 'jobsCosine/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				remote_job_folder = remote_base_folder + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

				if os.path.exists(local_job_folder+"timing.txt"): #Have we already copied this job over?
					print(f"Job already copied: {remote_job_folder}")
					continue

				#If you want to scp a whole folder recursivly it will copy the final folder 
				#which is stupid. To counter this, take out the last folder in local_path
				#since these are the same folder
				local_job_folder = '/'.join(local_job_folder.split('/')[:-2]) + '/'
				

				remote_file_exists = list_remote_files('Cosine', remote_job_folder+"timing.txt")
				if len(remote_file_exists) > 0 and "timing.txt" == remote_file_exists[0].split('/')[-1]:
					if not os.path.exists(local_job_folder): #folder doesnt exist locally so make the folder(s)
						os.makedirs(local_job_folder)
					print(f"Copying {remote_job_folder}")
					scp_transfer_with_ssh_config('Cosine',remote_job_folder,local_job_folder,recursive=True)


if __name__ == '__main__':
	main()