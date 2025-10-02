import numpy as np
import itertools

import matplotlib.pyplot as plt

import mdtraj as md

import time
import os
import sys

sys.path.insert(1, f'{os.getcwd()}/../utility')
import compute_observable_on_trjs_x01_v3


def rmsd_rmsf(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/aa_z_dist_2/"
    _savefilename = upperpath + savefilename

    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    #calculate protein center of mass (COM) by frame
    protein_atoms = trj.top.select("protein and name CA")
    prot_com = md.compute_center_of_mass(trj.atom_slice(protein_atoms))

    aa_inds = trj.top.select("resname ARAN and name C1")
    print(aa_inds)
    np.save(_savefilename+"-aaz" + ".npy", trj.xyz[:,aa_inds,2]-np.tile(prot_com[:,2], (16,1)).transpose())



compute_observable_on_trjs_x01_v3.compute_observable_x01(rmsd_rmsf)
