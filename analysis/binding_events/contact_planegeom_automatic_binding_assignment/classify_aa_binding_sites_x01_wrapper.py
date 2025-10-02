import numpy as np
import mdtraj as md
#from protein_ligand_contacts import get_contacts
#from protein_ligand_contacts import save_pdb_bfactors
import time
#from enspara.info_theory.exposons import exposons_from_sasas
import os
import sys

import classify_aa_binding_sites

sys.path.insert(1, f'{os.getcwd()}/../utility')

import load_all_trjs
#import load_anatale_trjs
import itertools
import count_lipids

import compute_observable_on_trjs_x01_v3


dummy_atom_queries = {"aac1":"resSeq 596", "ucp1":"resSeq 576"}

def closest_approaches(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/binding_planegeom"
    stride = 1

    #load whole trajectory once
    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")
    
    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices=ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    #set dummy atom coordinates to protein center of mass
    all_prot_indices = trj.top.select(f"protein and name CA")
    prot_com = md.compute_center_of_mass(trj.atom_slice(all_prot_indices))
    for t in range(trj.n_frames):
        for k in range(3):
            trj.xyz[t][-1][k] = prot_com[t][k]
            
    dummy_atom_query = dummy_atom_queries[protein]
    events_by_helix, events_by_helix_fa = classify_aa_binding_sites.binding_state(trj, dummy_atom_query, protein)

    #save contacts
    a_savefilename = f"{upperpath}/{savefilename}-plane-crossings"
    np.save(a_savefilename, events_by_helix.astype("uint8"))

    b_savefilename = f"{upperpath}/{savefilename}-plane-fa-crossings"
    np.save(b_savefilename, events_by_helix_fa.astype("uint8"))


compute_observable_on_trjs_x01_v3.compute_observable_x01(closest_approaches)
