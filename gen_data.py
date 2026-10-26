"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file goes through all folders in a specified job set and calculates the specified data for this run
for the folders which have completed jobs in them.

"""



import sys
import glob
import os
import re
import json
import csv
import numpy as np
from pathlib import Path


relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

# Path for the single global log file
GLOBAL_LOG = f"{project_path}../SpaceLab_data/logs/"

with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
path = input_json["data_directory"]

#from Suyama 2012
#Just geometric cross section, no normalization or nothing. 
def calc_geometric_cross_section(data_folder,size,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan

    if makeVisual:
        orientations = 1
    else:
        orientations = 30
    
    gcs = 0
    direction = (0,0,1)
    for _ in range(orientations):
        if not makeVisual:
            direction = u.make_rand_vector(3)
            gcs += u.calc_geometric_cross_section(pos,radius,direction=direction,write=makeVisual,directory=data_folder)
        else:
            gcs += u.calc_geometric_cross_section(pos,radius,direction=direction,mesh_factor=1.0,write=makeVisual,directory=data_folder)
    gcs /= orientations
    
    # normalization = np.pi*np.sum(np.power(radius,2))

    if makeVisual:
        R = np.sqrt(gcs/np.pi)
        COM = u.calcCOM(mass,pos)
        u.writeGCSRadius(R,COM,data_folder)

    return gcs


def calc_porosity_ch(data_folder,size,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan
    
    effective_radius_cubed = np.sum(np.power(radius,3))

    hull_vol = u.calc_hull_volume(pos,radius,write=makeVisual,directory=data_folder)

    #if the MVEE failed, return NAN
    if np.isnan(hull_vol):
        return np.nan

    porosity = 1-((4.0/3.0)*np.pi*effective_radius_cubed)/(hull_vol)

    return porosity

def calc_porosity_fee(data_folder, size, relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    # Load particle data
    pos, radius, mass, moi = u.get_data(data_folder, data_index=size, relax=relax)
    if pos is None:
        return np.nan

    # Compute reference volume
    effective_radius_cubed = np.sum(radius**3)

    # Compute ellipsoid + QC
    status, max_q, c, R, axes = u.calc_fully_enclosed_ellipsoid(pos, radius, write=makeVisual,directory=data_folder)

    # --- Logging ---
    logfile = GLOBAL_LOG+"ellipsoid_log.csv"
    log_exists = os.path.isfile(logfile)

    row = {
        "data_folder": data_folder,
        "size": size,
        "status": status,
        "max_q": max_q
    }

    with open(logfile, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row.keys())
        if not log_exists:
            writer.writeheader()
        writer.writerow(row)

    if status != "ok":
        return np.nan

    # Compute porosity
    porosity = 1 - effective_radius_cubed / (axes[0] * axes[1] * axes[2])

    return porosity


def calc_porosity_fes(data_folder,size,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan
        
    effective_radius_cubed = np.sum(np.power(radius,3))

    _,full_sphere_radius = u.calc_fully_enclosed_sphere(pos,radius,write=makeVisual,directory=data_folder)

    porosity = 1-effective_radius_cubed/np.power(full_sphere_radius,3)

    return porosity

def calc_asymmetry_parameter(data_folder,size,relax=False,**kwargs):
    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan
    # num_balls = data.shape[0]
    if "bigbox" in data_folder:
        boxsize = u.value_from_directory("B",data_folder)
        pos = u.unwrap_periodic_positions(pos,boxsize)

    vol = [(4*np.pi/3)*r**3 for r in radius]

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 
        
    principal_moi = u.get_principal_moi(mass,pos,radius)[::-1]
    alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

    # return alphai[0]/np.sqrt(alphai[1]*alphai[2]) #From Draine Sensitivity of Polarization to Grain Shape. I. Convex Shapes
    return np.sqrt(alphai[0]/(alphai[1]+alphai[2]-alphai[0])) #From Draine Sensitivity of Polarization to Grain Shape. II. Aggregates


def calc_stretch_parameter(data_folder,size,relax=False,**kwargs):
    data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if data is None:
        return np.nan
    # num_balls = data.shape[0]

    vol = [(4*np.pi/3)*r**3 for r in radius]

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

        
    principal_moi = u.get_principal_moi(mass,data,radius)[::-1]
    alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

    return alphai[1]/np.sqrt(alphai[0]*alphai[2])

def calc_porosity_abc(data_folder,size,relax=False, **kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if "bigbox" in data_folder:
        boxsize = u.value_from_directory("B",data_folder)
        pos = u.unwrap_periodic_positions(pos,boxsize)

    if pos is None:
        return np.nan

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

    a,b,c = u.calc_equivalent_ellipsoid_principal_axes(effective_radius,pos,mass,radius,data_folder,write=makeVisual)

    porosity = 1-(effective_radius**3/(a*b*c))

    return porosity

def calc_effective_radius(data_folder,size,relax=False, **kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

    return effective_radius

def calc_porosity_KBM(data_folder,size,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
    if pos is None:
        return np.nan

    effective_radius = np.power(np.sum(np.power(radius,3)),1/3)   

    RKBM = u.calc_gyration_radius(effective_radius,pos,mass,radius,data_folder,write=makeVisual)

    porosity = 1-np.power((effective_radius/RKBM),3)

    return porosity

def calc_number_of_contacts(data_folder,data_index=-1,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)

    return u.calc_max_number_of_contacts(data_folder,data_index,relax,makeVisual)

    



def calc_fractal_dimension(data_folder,size,relax=False,**kwargs):
    makeVisual = kwargs.get("makeVisual", False)
    print(
        f"Before aggregate {directory}: "
        f"max RSS = {u.max_rss_gib():.2f} GiB",
        flush=True,
    )
    overwrite_octree_data = True
    show_FD_plots = False
    o3dv = u.o3doctree(data_folder,ppb=20000,overwrite_data=overwrite_octree_data,index=size,Temp=-1,relax=relax)
    o3dv.make_tree()
    FD = o3dv.calc_fractal_dimension(show_graph=show_FD_plots)
    return FD



data_functions = [calc_porosity_abc,calc_porosity_KBM, \
                    calc_number_of_contacts,calc_fractal_dimension, \
                    calc_asymmetry_parameter,calc_stretch_parameter, \
                    calc_porosity_fee,calc_porosity_fes, \
                    calc_porosity_ch,calc_geometric_cross_section, \
                    calc_effective_radius]
data_headers = [i.__name__[5:] for i in data_functions]

def calc_from_size(size,directory,existing_headers,existing_values,requested_headers,relax=False,overwrite=False,makeVisual=False):
    headers = []
    values = []

    for h,header in enumerate(data_headers):
        if not overwrite:
            #we already have this header include it
            if header in existing_headers:
                index = existing_headers.index(header)
                headers.append(header)
                #if the value is nan then we need to recalculate
                if existing_values[index] == "nan":
                    existing_values[index] = str(data_functions[h](directory,size,relax,makeVisual=makeVisual))
                values.append(existing_values[index])
            #we dont have this header and we want it, so calculate
            elif (header not in existing_headers and header in requested_headers):
                headers.append(header)
                val = np.nan
                for _ in range(10):
                    val = data_functions[h](directory,size,relax,makeVisual=makeVisual)
                    if not np.isnan(val):
                        break
                else:
                    print(f"Non nan value not calculated for directory: {directory}")
                values.append(str(val))
        else:
            #we already have this header and didn't ask for it, so just include it
            if header in existing_headers and header not in requested_headers:
                index = existing_headers.index(header)
                headers.append(header)
                values.append(existing_values[index])
            #if overwrite, we want to calculate reguardless, if header is requested
            elif header in requested_headers:
                headers.append(header)
                val = np.nan
                for _ in range(10):
                    val = data_functions[h](directory,size,relax,makeVisual=makeVisual)
                    if not np.isnan(val):
                        break
                else:
                    print(f"Non nan value not calculated for directory: {directory}")
                values.append(str(val))

    return headers,values

def extract_projectiles(filename):
    """
    Find projectile index/folder pairs for simulations that reached
    'Simulation complete!' before another projectile was started.

    Returns
    -------
    completed : list of tuples
        [(index, folder), ...]
    """

    pattern = re.compile(
        r"Getting projectile of index\s+(\d+)\s+from:\s+(.+)"
    )

    completed = []
    seen = set()
    pending = None

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            match = pattern.search(line)

            if match:
                index = int(match.group(1))
                folder = match.group(2).strip()

                # # Keep only the path after SpaceLab_data/
                # folder = folder.split("SpaceLab_data/", 1)[-1]
                # Remove double //
                folder = folder.replace("//","/")

                # Replace any previous unfinished projectile
                pending = (index, folder)
                continue

            if "Simulation complete!" in line and pending is not None:
                # Only add this exact (index, folder) pair once
                if pending not in seen:
                    completed.append(pending)
                    seen.add(pending)

                pending = None

    return completed

def find_projectiles():
    with open(project_path+"default_files/default_input.json",'r') as fp:
        input_json = json.load(fp)
    
    path = input_json["data_directory"]

    simfiles = ["sim_err.log","sim_errors.txt"]

    all_completed = []
    seen = set()

    # #CAS and wCAS (BAPA)
    # A = [i for i in range(30)]
    # M = [3,5,10,15,20,30,50,60,75,150]
    # N = [300]
    # T = [1000]

    # for a in A:
    #     for m in M:
    #         for n in N:
    #             for t in T:
    #                 simpath = path + f"jobs/BAPA_{a}/M_{m}/N_{n}/T_{t}/"
    #                 if (os.path.exists(simpath + simfiles[0])):
    #                     results = extract_projectiles(simpath + simfiles[0])
    #                 elif (os.path.exists(simpath + simfiles[1])):
    #                     results = extract_projectiles(simpath + simfiles[1])
    #                 else:
    #                     print(f"Sim output does not exist for {simpath}.")
    #                     exit(0)
    #                 for pair in results:
    #                     if pair not in seen:
    #                         all_completed.append(pair)
    #                         seen.add(pair)

    # #CNP (CBAPA))
    # A = [i for i in range(30)]
    # M = [3,5,10,15,20,30,50]
    # C = 30
    # N = [m*C for m in M]
    # T = [1000]

    # for a in A:
    #     for i,m in enumerate(M):
    #         n = N[i]
    #         for t in T:
    #             simpath = path + f"jobs/CBAPA_{a}/M_{m}/N_{n}/T_{t}/"
    #             if (os.path.exists(simpath + simfiles[0])):
    #                 results = extract_projectiles(simpath + simfiles[0])
    #             elif (os.path.exists(simpath + simfiles[1])):
    #                 results = extract_projectiles(simpath + simfiles[1])
    #             else:
    #                 print(f"Sim output does not exist for {simpath}.")
    #                 exit(0)
    #             for pair in results:
    #                 if pair not in seen:
    #                     all_completed.append(pair)
    #                     seen.add(pair)

    #CFS (DBAPA))
    A = [i for i in range(30)]
    M = 20
    C = [1,2,3,4,5,6,10,15,20,30]
    N = [c*M for c in C]
    T = [1000]

    for a in A:
        for i,n in enumerate(N):
            m = M
            for t in T:
                simpath = path + f"jobs/DBAPA_{a}/M_{m}/N_{n}/T_{t}/"
                if (os.path.exists(simpath + simfiles[0])):
                    results = extract_projectiles(simpath + simfiles[0])
                elif (os.path.exists(simpath + simfiles[1])):
                    results = extract_projectiles(simpath + simfiles[1])
                else:
                    print(f"Sim output does not exist for {simpath}.")
                    exit(0)
                for pair in results:
                    if pair not in seen:
                        all_completed.append(pair)
                        seen.add(pair)
    
    all_completed.sort(key=lambda x: (x[0], x[1]))
    for comp in all_completed:
        print(comp)

def main():


    # find_projectiles()
    # exit(0)

    # print(calc_porosity_abc("/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_0/N_300/T_1000/",300,True))
    # with open(project_path+"default_files/default_input.json",'r') as fp:
    #     input_json = json.load(fp)
    
    # path = input_json["data_directory"]



    # data_file = "job_data.csv" 
    data_file = "job_data_temp.csv" 



    N = [1000]
    N = [300] 
    C = 30
    M = [1,3,5,10,15,20,30,50,60,75,100,150]
    # M = [3,15,100]

    #list of the functions that calculate the data you want
    #if adding to this list, name your function calc_*header_name*
    #where *header_name* is what you want this data to be called
    #in the header of the job_data.csv file. The function should 
    #only take in the directory the data is in and the size of which 
    #to calculate as an input.
    #It should return a single data value.
    # bool_headers = [1,1,0,0,1,0,0,0]

    bool_headers = [1,1,0,0,1,0,0,0,0,0,0]
    bool_headers = [1,0,0,0,1,0,0,0,0,0,0]
    bool_headers = [1,1,0,0,1,1,0,0,0,1,0]

    # requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
    requested_data_headers = [data_headers[i] for i in range(len(data_headers)) if bool_headers[i]]

    overwritedata = False
    makeVisual = False

    print("=========================START GEN DATA=========================")

    for n_i,n in enumerate(N):
    
        data_folders = []
        # data_folders += [path + 'jobs/DBAPA_*']
        # data_folders += [path + 'jobs/CBAPA_*']
        data_folders += [path + 'jobs/BAPAWELD_*']
        # data_folders += [path + 'jobs/BAPA_*']
        # data_folders += [path + 'jobsCosine/lognorm_*/N_300/T_1000/']
        # data_folders += ["/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_branch/SpaceLab_data/" + 'jobs/bigbox_*']


        # #BAPA
        # data_folders += [path + 'jobsCosine/lognorm_*/N_30/T_1000/']
        # data_indices = [3,5,10,15]
        # data_folders += [path + 'jobsCosine/lognorm_*/N_300/T_1000/']
        # data_indices = [3,5,10,15,20,30,50,60,75,150]
        # #CBAPA
        # data_folders += [path + 'jobsCosine/lognorm_*/N_300/T_1000/']
        # data_indices = [3,5,10,15,20,30,50]
        # #DBAPA
        # data_folders += [path + 'jobsCosine/lognorm_*/N_300/T_1000/']
        # data_indices = [20]


        possible_dirs = []
        for data_folder in data_folders:
            # print(data_folder)
            possible_dirs.extend(u.get_directores_containing(data_folder,["timing.txt"]))
        # print(possible_dirs)

        #list of intermediate sizes to calculate data for.

        # requested_sizes = list(range(30,301))
        # requested_sizes = [[u.find_max_index(directory)] for directory in possible_dirs]
        # requested_sizes = [[n] for directory in possible_dirs]

        for d_i,directory in enumerate(possible_dirs):
            CBAPA = ("CBAPA" in directory) #this is for if we want to use CBAPA since the N and M are different than the usual constant N 
            DBAPA = ("DBAPA" in directory) #this is for if we want to use DBAPA since the N and M are different than the usual constant N 
            #en is the n value to calculate data for
            if CBAPA:
                # n = u.find_max_index(directory)
                m = u.value_from_directory("M",directory)
                en = C*m
            elif DBAPA:
                m = u.value_from_directory("M",directory)
                en = u.value_from_directory("N",directory)
            else:
                en = n


            requested_sizes = [en]
            # requested_sizes = data_indices
            relax = ("relax" in directory)
            # relax = False
            # print(f"relax: {relax}")
            rel = ""
            if relax:
                rel = "relax_"
            # directory = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsNovus/constrelax_10/N_30/T_1000/"
            print(f"Generating data for: {directory}{rel}{data_file}")

            if os.path.exists(directory+"timing.txt"):
                existing_data = []
                if os.path.exists(directory+f"{rel}{data_file}"):
                    with open(directory+f"{rel}{data_file}",'r') as fp:
                        existing_data = fp.readlines()
                existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]

                #This will be all the data we want to write in the end
                lines = []
                
                #contains both the sizes we want and the sizes we already have
                all_sizes = sorted(set(existing_sizes+requested_sizes))
                # all_sizes = sorted(set(existing_sizes+requested_sizes[d_i]))
                # all_sizes = [300]
                # print(all_sizes)
                for size in all_sizes:

                    existing_headers_for_size = []
                    existing_values_for_size = []
                    if size in existing_sizes:
                        #calculate the index of the existing size
                        index = existing_sizes.index(size)*4
                        existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
                        existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")

                    #If this size wasn't requested then we don't want to overwrite it
                    #Or if overwrite is true then overwrite anyway
                    if overwritedata or size in requested_sizes:
                    # if overwritedata or size in requested_sizes[d_i]:
                        headers,values = calc_from_size(size,directory,existing_headers_for_size,
                                                        existing_values_for_size,requested_data_headers,
                                                        relax,overwritedata,makeVisual)
                        print("===========================================")
                        # print(size,directory)
                        print(headers)
                        print(values)
                    else:
                        headers = existing_headers_for_size
                        values = existing_values_for_size


                    lines.append(f"N={size}\n")
                    lines.append(','.join(headers)+'\n')
                    lines.append(','.join(values)+'\n')
                    lines.append('\n')


                with open(directory+f"{rel}{data_file}",'w') as fp:
                    fp.writelines(lines)
                # print(lines)

            else:
                print(f"Job is not finished in {directory}")
            

def add_fractal_dimension(
    folder,
    source_file="job_data.csv",
    target_file="job_data_temp.csv",
):
    folder = Path(folder)
    source = folder / source_file
    target = folder / target_file

    if not source.is_file() or not target.is_file():
        print(f"Missing file in {folder}")
        return 0

    def read_blocks(filename):
        blocks = {}

        lines = filename.read_text().splitlines()

        for i, line in enumerate(lines):
            if not line.startswith("N="):
                continue

            n = int(line.split("=", 1)[1])
            headers = lines[i + 1].split(",")
            values = lines[i + 2].split(",")

            blocks[n] = dict(zip(headers, values))

        return blocks

    source_blocks = read_blocks(source)

    lines = target.read_text().splitlines()
    updated = 0

    for i, line in enumerate(lines):
        if not line.startswith("N="):
            continue

        n = int(line.split("=", 1)[1])

        if n not in source_blocks:
            continue

        if "fractal_dimension" not in source_blocks[n]:
            continue

        headers = lines[i + 1].split(",")
        values = lines[i + 2].split(",")

        fd = source_blocks[n]["fractal_dimension"]

        if "fractal_dimension" in headers:
            index = headers.index("fractal_dimension")
            values[index] = fd
        else:
            headers.append("fractal_dimension")
            values.append(fd)

        lines[i + 1] = ",".join(headers)
        lines[i + 2] = ",".join(values)

        updated += 1

    if updated:
        target.write_text("\n".join(lines) + "\n")

    return updated

def add_fractal_dimension_to_all(root):
    root = Path(root)

    total_files = 0
    total_blocks = 0

    for target in root.rglob("job_data_temp.csv"):
        updated = add_fractal_dimension(target.parent)

        if updated:
            total_files += 1
            total_blocks += updated

    print(f"Updated {total_blocks} blocks in {total_files} files.")


if __name__ == '__main__':

    # for run in ["BAPA","CBAPA","BAPAWELD","DBAPA"]:
    #     for folder in Path(Path(path) / Path("jobs")).glob(f"{run}_*"):
    #         add_fractal_dimension_to_all(folder)

    main()

