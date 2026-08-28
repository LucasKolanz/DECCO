"""
This file was originally written for SpaceLab/DECCO to check the seeds of all jobs in folders matching a specified pattern.

Author: Lucas Kolanz

This file is meant to go through all jobs contained in a specified file pattern and print out their seeds. All of these seeds should
be different, and this file is meant to help verify that.

TODO: use the seeds from input.json, not seedFile.txt. Save all seeds in a list to verify they are all unique.

"""





import os
import json
import glob
import numpy as np
import subprocess
import sys
import traceback
from pathlib import Path

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u



def calc_agg_radii(pos,COM):
	
	radii = 0
	for p in pos:
		r = np.linalg.norm(p-COM)
		
		if r > radii:
			radii = r 

	return radii



def main():


	# Open SpaceLab default file for directory information
	with open(project_path + "default_files/default_input.json", "r") as fp:
		input_json = json.load(fp)

	data_directory = Path(input_json["data_directory"])

	jobs_root = Path(input_json["data_directory"]) / "jobs"

	job_globs = [
		"BAPA_*/M_*/N_*/T_*",
		"CBAPA_*/M_*/N_*/T_*",
		"DBAPA_*/M_*/N_*/T_*",
		"BAPAWELD_*/M_*/N_*/T_*",
	]

	folders = sorted(
		folder
		for pattern in job_globs
		for folder in jobs_root.glob(pattern)
		if folder.is_dir() and (folder / "timing.txt").is_file()
	)

	print(f"Found {len(folders)} completed folders")

	aggradii = []
	monomerradii = []
	agg_folders = []

	for folder in folders:
		agg_folders.append(str(folder) + "/")
		try:
			pos, radii, mass, _ = u.get_data(str(folder) + "/", relax=False)
			# print(pos.shape)
			aggradii.append(calc_agg_radii(pos, u.calcCOM(mass, pos)))

			monomerradii.append(max(radii))
			monomerradii.append(min(radii))

		except Exception:
			print(f"Error processing: {folder}")
			traceback.print_exc()
			exit(0)

	aggradii = np.asarray(aggradii, dtype=np.float64)

	if aggradii.size > 0:
		print(f"Processed {len(aggradii)} folders")
		print(f"Maximum aggregate radius: {np.nanmax(aggradii)} cm")
		print(f"Maximum aggregate folder: {agg_folders[np.nanargmax(aggradii)]} cm")
		print(f"Minimum aggregate radius: {np.nanmin(aggradii)} cm")
		print(f"Minimum aggregate folder: {agg_folders[np.nanargmin(aggradii)]} cm")

		print(f"Maximum monomer radius: {np.nanmax(monomerradii)} cm")
		print(f"Minimum monomer radius: {np.nanmin(monomerradii)} cm")
	else:
		print("No matching folders were successfully processed.")

if __name__ == '__main__':
	main()