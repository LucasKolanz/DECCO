Welcome to SpaceLab!

SpaceLab is a discrete element code for simulating cosmic collisions. These collisions can be at the scale of micrometers for simulating dust to the scale of kilometers for simulating asteroid collisions. 

Dependencies:
	hdf5
	nhloman json


No installation is needed for SpaceLab. All you need is a g++ compiler and you should be good to go.

To get started, all you need to do is run one of the "make_*.py" files and it will compile and run the specified jobs. These "make_*.py" files work by creating an input.json file for the simulation inputs, copying necessary files to the correct directores, and then compiling and running the individual jobs. Depending on your system, you may need to edit the Makefile provided.

To initialize the submodule in the git folder run the follow command
	python3 config.py
It will configure the external git repo and set some default project variables

File/Directory scheme:
When you make a new job set with a make_*.py file, there are several ways of organizing the resulting data. The first way to distinguish specific jobs is by setting the "job_set_name" variable in the make_*.py file. This will create a new folder in the form:
	/*SpaceLabDir*/jobs/job_set_name
This is useful both for small tests, where the output can go directly into the job_set_name folder, as well as large sets of jobs that can then have a further folder hierarchy to distinguish between specific jobs. In this case the make_*.py file may need to be edited to suit your job set's specific needs. For example, if I want a job set that has simulations for a set of temperatures for various aggregate final sizes, I might use the folder scheme:
	/*SpaceLabDir*/jobs/job_set_name/N_{}/T_{}/
Currently, if multiple runs of the same conditions (in this example, the conditions of N and T) then multiple job_set_name folders are made, distinguished by an index immediatly following job_set_name. For example:
	/*SpaceLabDir*/jobs/job_set_name0/N_{}/T_{}/
	/*SpaceLabDir*/jobs/job_set_name1/N_{}/T_{}/
There is no reason that this couldn't be changed to: 
	/*SpaceLabDir*/jobs/job_set_name/index/N_{}/T_{}/
Where index is the distinguishing feature.




Data files:
Data files in DECCO are in either csv or hdf5 format. There are three types of data file in DECCO
	constants 
		rows represent different balls
		columns are radius, mass, moment of inertia
	simData 
		rows represent timesteps
		columns are posx[0],posy[0],posz[0],wx[0],wy[0],wz[0],wmag[0],velx[0],vely[0],velz[0],bound[0],posx[1],posy[1],...
		(where 0,1,2,... are the different balls)
	energy
		rows represent timesteps
		columns are time, PE, KE, Etot, p, L