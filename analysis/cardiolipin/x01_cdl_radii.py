import mdtraj as md
import numpy as np
import os
import sys
import itertools

sys.path.insert(1, f'{os.getcwd()}/../utility')
import compute_observable_on_trjs_x01_v3


def calculate_radial_distances(trj, protein):

    #---------------------------part 1: trajectory preparation-------------------------------

    #calculate protein center of mass (COM) by frame
    protein_atoms = trj.top.select("protein and name CA")
    prot_com = md.compute_center_of_mass(trj.atom_slice(protein_atoms))

    #copy the trajectory (there were already no hydrogens) prior to editing
    trj_xy = trj.atom_slice(trj.top.select("not element H")) 

    #select atom to place at protein COM (a K+/Cl- ion in this case)
    dummy_atom_query = {"aac1":"resSeq 590", "ucp1":"resSeq 576"}
    dummy_ind = trj_xy.top.select(dummy_atom_query[protein])
    
    #set the coordinates of an atom with no other use in this analysis 
    # to equal the protein center of mass (COM) so that 
    # mdtraj minimum image distance calculation can be applied
    for t in range(trj_xy.n_frames):
        for k in range(3):
            trj_xy.xyz[t][dummy_ind[0]][k] = prot_com[t][k]
    
    #project everything to the z=0 plane,
    # which must be parallel to the average plane of the membrane since the box is periodic
    trj_xy.xyz[:,:,2] = 0 


    #---------------------------part 2: calculate protein-lipid distances-------------------------------

    #this could be any molecule which lives in the plane of the membrane
    lipid_atoms = trj_xy.top.select("resname TLCL2 and name C2")

    #calculate pairwise distances
    # doing it this way (i.e. using MDtraj's minimum image distance calculation for the phosphates) 
    # is important because reimaging centers molecules so the phosphate may be over the periodic boundary
    index_pairs = list(itertools.product(dummy_ind, lipid_atoms))
    lipid_radii = md.compute_distances(trj_xy, index_pairs, periodic=True, opt=True).reshape((trj_xy.n_frames, len(lipid_atoms))) #probably needlessly slow

    return lipid_radii


def x01_cdl_radii(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/cdl"

    #load whole trajectory once
    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")
    
    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices=ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    lipid_radii = calculate_radial_distances(trj, protein)

    np.save(f"{upperpath}/{savefilename}-cdlradii.npy", lipid_radii)


compute_observable_on_trjs_x01_v3.compute_observable_x01(x01_cdl_radii)
