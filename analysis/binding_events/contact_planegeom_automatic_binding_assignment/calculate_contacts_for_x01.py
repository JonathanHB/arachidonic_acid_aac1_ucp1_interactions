import numpy as np
import mdtraj as md
#import protein_ligand_contacts

#import load_trjs_compute_dists_concurrent
# from protein_ligand_contacts import get_trj_contacts
# from protein_ligand_contacts import save_pdb_bfactors
import time
#from enspara.info_theory.exposons import exposons_from_sasas

import os
import sys

import matplotlib.pyplot as plt

sys.path.insert(1, f'{os.getcwd()}/../utility')
#import importlib

import compute_observable_on_trjs_x01_v3

#import compute_observable_on_trjs
#import generate_query
import count_lipids
import itertools


#----------------------------------------------------------------------------------------------------
#                               residues for query construction

c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}

terminal_residues = {
    "aac1":
    [[9,   36],
    [74,  94],
    [114, 141],
    [177, 197],
    [211, 238],
    [274, 294]],
    "ucp1":
    [[14,  41],
    [78,  98],
    [114, 141],
    [177, 197],
    [213, 240],
    [271, 291]]
}


#----------------------------------------------------------------------------------------------------
#                                       protein queries

def resqueries(protein):

    helix_resseqs_all = []
    resqueries = []
    rq_all = ""
    
    for rt in terminal_residues[protein]:
        helix_resseqs = [i for i in range(rt[0], rt[1]+1)]
        helix_resseqs_all += helix_resseqs
        rq_all += "+".join([str(i) for i in helix_resseqs])
        #print("color red, resi " + "+".join([str(i) for i in helix_resseqs]))
        resqueries.append(" or ".join([f"resSeq {i}" for i in helix_resseqs]))

    #print(rq_all)
    
    return resqueries, helix_resseqs_all


#----------------------------------------------------------------------------------------------------
#                                       lipid/fa queries

def get_ligand_query(protein, residue, leaflet, refpath):
    #parameters
    #names = {"ARAN":["C1", "C5", "C10", "C15", "C20"], "POPC":["C21", "C31"]}
    names = {"ARAN":["O1","O2"] + ["C"+str(ci) for ci in range(1,21)], "POPC":["C21", "C31"]}
    namequery = " or ".join(["name "+ name for name in names[residue]])

    #select lipids or fatty acids
    if residue == "POPC":
        init_rsqs = {"aac1":1, "ucp1":9}
        upper_leaflet_ref_query = f"resSeq {init_rsqs[protein]} and name CA"
        upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet = count_lipids.count_lipids(refpath, upper_leaflet_ref_query, residue, "P")
        if leaflet == "upper":
            leaflet_rsq = upperleaflet_rsq
            nleaflet = upperleaflet
        elif leaflet == "lower":
            leaflet_rsq = lowerleaflet_rsq
            nleaflet = lowerleaflet
        else: 
            raise ValueError("leaflet must be upper or lower")
        
        c_l_query = " or ".join(["resSeq "+str(i) for i in leaflet_rsq])
        ligand_query = f"resname {residue} and {namequery} and ({c_l_query})"

    elif residue == "ARAN":
        c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
        invstring = ""
        if leaflet == "lower":
            invstring = " not "
        c_aran_query = " or ".join(["resSeq "+str(i) for i in c_aran[protein]])
        ligand_query = f"resname {residue} and ({namequery}) and {invstring}({c_aran_query})"
        nleaflet = 8 #number of ARAN in each leaflet

    else:
        raise ValueError("residue must be POPC or ARAN")

    return names[residue], f"({ligand_query})", nleaflet


#----------------------------------------------------------------------------------------------------
#                                       contact calculation

def get_trj_contacts(trj, protein_query, ligand_query, contact_threshold):
    
    #select protein and ligand atoms
    indices1 = trj.top.select(ligand_query)
    indices2 = trj.top.select(protein_query)

    # list of atom pairs
    index_pairs = list(itertools.product(indices1, indices2))
    
    # calculate pairwise distances
    # periodic = True is an issue for rcsb pdb structures with unit cell information if the residues are far enough apart
    distances = md.compute_distances(trj, index_pairs, periodic=True, opt=True).reshape(trj.n_frames, len(indices1), len(indices2))

    #calculate the minimum distance from each protein atom to any ligand atom
    mindists = np.min(distances, axis=1)
    contacts = np.where(mindists < contact_threshold, 1, 0)
    contact_frequency = np.mean(contacts, axis=0)

    return indices2, contacts, contact_frequency


#----------------------------------------------------------------------------------------------------
#                                       contact calculation

def save_contacts(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/contacts"

    #load whole trajectory once
    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    #get protein query
    rq, hr = resqueries(protein)
    helixresis_all = "or".join(rq)
    protein_query = f"protein and ({helixresis_all}) and not element H"

    #loop over relevant molecules and leaflets
    #for ligand in ["ARAN", "POPC"]:
    for ligand in ["ARAN"]:
        for leaflet in ["upper", "lower"]:

            #get lipid/fa query
            name, ligand_query, nlipids = get_ligand_query(protein, ligand, leaflet, top_path)

            #break up trajectories too large to fit in memory (on x01)
            #step = 5000
            #for i in range(0, trj.n_frames, step):

            #endframe = min(i+step, trj.n_frames)
            #print(f"processing frames {i} to {endframe}")

            #calculate contacts
            contact_dist_threshold = 0.7
            indices, contacts, contact_freq = get_trj_contacts(trj, protein_query, ligand_query, contact_dist_threshold)

            #save contacts
            freq_savefilename = f"{upperpath}/{savefilename}-{leaflet}-{ligand}-contact-freq"
            np.save(freq_savefilename, contact_freq)
        
            t_savefilename = f"{upperpath}/{savefilename}-{leaflet}-{ligand}-contacts"
            np.save(t_savefilename, contacts)



def save_contacts_by_fa(trj_path, top_path, protein, savefilename):

    contact_dist_threshold = 0.5
    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/contacts"

    #load whole trajectory once
    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    #get protein query
    rq, hr = resqueries(protein)
    helixresis_all = "or".join(rq)
    protein_query = f"protein and ({helixresis_all}) and not element H"

    #loop over relevant molecules and leaflets
    #for ligand in ["ARAN", "POPC"]:
    for ligand in ["ARAN"]:
        for leaflet in ["upper"]:

            contacts_by_aa = []
            #get lipid/fa query
            #name, ligand_query, nlipids = get_ligand_query(protein, ligand, leaflet, top_path)
            c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
            for aa_resseq in c_aran[protein]:
                #invstring = ""
                # if leaflet == "lower":
                #     invstring = " not "
                #c_aran_query = " or ".join(["resSeq "+str(i) for i in c_aran[protein]])
                ligand_query = f"resname {ligand} and resSeq {aa_resseq} and not element H"

                #calculate contacts
                indices, contacts, contact_freq = get_trj_contacts(trj, protein_query, ligand_query, contact_dist_threshold)
                contacts_by_aa.append(contacts)

            contacts_by_aa = np.stack(contacts_by_aa).astype(np.uint8)

            #save contacts
            savefilename = f"{upperpath}/{savefilename}-{leaflet}-{ligand}-contacts-byaa"
            np.save(savefilename, contacts_by_aa)
        
            # t_savefilename = f"{upperpath}/{savefilename}-{leaflet}-{ligand}-contacts"
            # np.save(t_savefilename, contacts)



compute_observable_on_trjs_x01_v3.compute_observable_x01(save_contacts_by_fa)