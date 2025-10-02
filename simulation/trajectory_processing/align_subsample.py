import os
import numpy as np

import sys
sys.path.insert(1, f'{os.getcwd()}/../../analysis/utility')
import compute_observable_on_trjs_x01_v3 

def subsample_frames(trj_path, top_path, protein, savefilename):

    outpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories_s10/"
    tprpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories/wynton/"

    g1 = os.system(f"echo 1 1 0 | gmx trjconv -f {trj_path} -s {tprpath+protein}/input/{protein}_500ns.tpr -o {outpath+savefilename}-aligned-s20.xtc -fit rot+trans -center -skip 20")


    # waterpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/water_wires_2/" + savefilename
    # #waterpath = waterupperpath + savefilename

    # outpath_u = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories_wire_s10/"
    # outpath = outpath_u + savefilename

    # wire_gap = np.load(f"{waterpath}-gaps.npy")
    # #print(wire_gap.shape)
    # frame_inds = np.argwhere(wire_gap < 0.33).flatten()
    # frame_ind_string = ' '.join([str(f+1) for f in frame_inds])

    # #print(frame_inds)

    # framefn = "/".join([outpath_u, savefilename.split("/")[0], savefilename.split("/")[1], savefilename.split("/")[2], f"{savefilename.split('/')[3][-5:]}-index.ndx"])
    # if os.path.exists(framefn):
    #     os.remove(framefn)

    # #https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2010-July/052593.html
    # os.system(f"echo [ frames ] >> {framefn}")
    # os.system(f"echo {frame_ind_string} >> {framefn}")

    # outpath_u = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories_wire_s10/"
    # outpath = outpath_u + savefilename

    # if len(frame_inds) > 0:
    #     os.system(f"echo 0 | gmx trjconv -f {trj_path} -s {top_path} -fr {framefn} -o {outpath}-wireframes.xtc")
    #     os.system(f"echo 0 | gmx trjconv -f {outpath}-wireframes.xtc -skip 10 -o {outpath}-wireframes-s10.xtc")


compute_observable_on_trjs_x01_v3.compute_observable_x01(subsample_frames)
