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

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/rmsd_rmsf/"
    _savefilename = upperpath + savefilename

    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    rmsd_atom_query = "protein and name CA"

    #also used for rmsf
    protein_to_inputfolder = {"aac1":"aac1-eq-run20-input", "ucp1":"ucp1-eq-run21-input"}
    rmsdref = md.load(f"/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/initial_structures/{protein_to_inputfolder[protein]}/step5_input.gro")
    rmsdref = rmsdref.atom_slice(rmsdref.top.select(rmsd_atom_query))

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select(rmsd_atom_query)
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    trj.center_coordinates()

    rmsd = md.rmsd(trj, rmsdref, 0, precentered=True)
    np.save(_savefilename+"-rmsd" + ".npy", rmsd)

    rmsf = md.rmsf(trj, rmsdref, 0, precentered=True)
    np.save(_savefilename+"-rmsf" + ".npy", rmsf)



compute_observable_on_trjs_x01_v3.compute_observable_x01(rmsd_rmsf)
