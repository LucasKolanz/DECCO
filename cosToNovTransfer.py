from paramiko import SSHClient, SSHConfig, AutoAddPolicy
from scp import SCPClient, SCPException
import os
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def ensure_remote_directory_exists(hostname, directory):
	"""
	Ensures that a directory exists on the remote server. If the directory does not exist, it is created.
	
	:param ssh: An established SSHClient connection.
	:param directory: The path of the directory to check/create on the remote server.
	"""

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
	
	command = f"mkdir -p {directory}"
	stdin, stdout, stderr = ssh.exec_command(command)
	exit_status = stdout.channel.recv_exit_status()  # Blocking call
	if exit_status == 0:
		print(f"Directory ensured: {directory}")
	else:
		print(f"Error ensuring directory {directory}. Exit status: {exit_status}")

	ssh.close()

#IF files is empty then just copy remote_file_path, otherwise loop through and copy all files from remote_file_path
def scp_transfer_with_ssh_config(hostname, source_file_path, dest_path, scp_command ,recursive=False): 
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
			# print(f"Attempting to copy {source_file_path} to {dest_path}")
			# print(recursive)
			# print(source_file_path)
			# print(dest_path)
			if scp_command == "get":
				scp.get(source_file_path, dest_path,recursive=recursive)
			elif scp_command == "put":
				scp.put(source_file_path, dest_path,recursive=recursive)
			else:
				print(f"ERROR: scp command '{scp_command}' doesn't exist. Should be either 'get' for remote to local or 'put' for local to remote.")

			file_found = True
		except SCPException as e:
			# Check if the exception is due to a missing file
			print(e)
			
			


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


	cosine_data_directory = '/home/physics/kolanzl/SpaceLab_data/'
	novus_data_directory = '/home/kolanzl/novus/kolanzl/SpaceLab_data/'

	#load default input file
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	local_data_directory = input_json["data_directory"]


	job_set_name = "lognorm"

	attempts = [i for i in range(30)]
	# attempts = [0,29]
	attempts_300 = attempts


	N = [30,100,300]
	# N=[300]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	for n in N:
		for Temp in Temps:
			temp_attempt = attempts
			if n == 300:
				temp_attempt = attempts_300
			for attempt in temp_attempt:
				local_job_directory = local_data_directory + 'tempJobs/'+job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				novus_job_directory = novus_data_directory  + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'
				cosine_job_directory = cosine_data_directory + 'jobs/' + job_set_name + str(attempt) + '/'\
							+ 'N_' + str(n) + '/' + 'T_' + str(Temp) + '/'

				# if os.path.exists(novus_job_directory+"timing.txt"): #Have we already copied this job over?
				# 	print(f"Job already copied: {novus_job_directory}")
				# 	continue

				#If you want to scp a whole folder recursivly it will copy the final folder 
				#which is stupid. To counter this, take out the last folder in local_path
				#since these are the same folder
				short_cosine_job_directory = '/'.join(cosine_job_directory.split('/')[:-2]) + '/'
				short_novus_job_directory = '/'.join(novus_job_directory.split('/')[:-2]) + '/'
				short_local_job_directory = '/'.join(local_job_directory.split('/')[:-2]) + '/'
				

				cosine_file_exists = list_remote_files('Cosine', cosine_job_directory+"timing.txt")
				if len(cosine_file_exists) > 0 and "timing.txt" == cosine_file_exists[0].split('/')[-1]: #if timing is there the sim is finished and doesn't need copying over. So only copy timing.txt so we know that one is done
					cosine_job_directory += "timing.txt"
					short_local_job_directory = local_job_directory
					local_job_directory = local_job_directory+"timing.txt"
					# print(f"cosine_job_directory: {cosine_job_directory}")
					# print(f"short_local_job_directory: {short_local_job_directory}")
					# print(f"local_job_directory: {local_job_directory}")
				
				if not os.path.exists(local_job_directory): #folder doesnt exist locally so make the folder(s)
					os.makedirs(local_job_directory.strip("timing.txt"))
				print(f"Copying {cosine_job_directory} to {novus_job_directory} via {local_job_directory}")
				scp_transfer_with_ssh_config('Cosine',cosine_job_directory,short_local_job_directory,scp_command='get',recursive=True)
				ensure_remote_directory_exists('Novus',novus_job_directory)
				scp_transfer_with_ssh_config('Novus',local_job_directory,novus_job_directory,scp_command='put',recursive=True)


if __name__ == '__main__':
	main()