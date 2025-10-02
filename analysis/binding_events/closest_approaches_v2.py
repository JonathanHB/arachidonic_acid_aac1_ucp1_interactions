import numpy as np
import mdtraj as md
#from protein_ligand_contacts import get_contacts
#from protein_ligand_contacts import save_pdb_bfactors
import time
#from enspara.info_theory.exposons import exposons_from_sasas
import os
import sys

sys.path.insert(1, f'{os.getcwd()}/../utility')

import load_all_trjs
#import load_anatale_trjs
import itertools
import count_lipids

import compute_observable_on_trjs_x01_v2


def get_distance_matrix(trj, protein_indices, ligand_indices, periodic=True):

    #select protein and ligand atoms
    indices1 = protein_indices
    indices2 = ligand_indices

    # list of atom pairs
    index_pairs = list(itertools.product(indices1, indices2))

    # calculate pairwise distances
    # periodic = True is an issue for rcsb pdb structures with unit cell information if the residues are far enough apart
    distances = md.compute_distances(trj, index_pairs, periodic=periodic, opt=True).reshape(trj.n_frames, len(indices1), len(indices2))

    #minimum distance from any protein atom to any ligand atom in each frame
    return np.min(distances, axis = (1,2))


def resqueries(protein):

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
    
    resqueries = []
    for rt in terminal_residues[protein]:
        helix_resseqs = [i for i in range(rt[0], rt[1]+1)]
        #print("color red, resi " + "+".join([str(i) for i in helix_resseqs]))
        resqueries.append(" or ".join([f"resSeq {i}" for i in helix_resseqs]))
    
    return resqueries


def get_residue_info(trj, protein, refpath, residue, leaflet):

    charge_center_ind = {"ARAN":2, "POPC":6} #index of atom within molecule
    ref_atoms = {"ARAN":"C1", "POPC":"P"} #atoms to use for head groups
    init_rsqs = {"aac1":1, "ucp1":9} #residues on C side of protein

    #identify residue in the desired leaflet
    upper_leaflet_ref_query = f"resSeq {init_rsqs[protein]} and name CA"
    upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet = count_lipids.count_lipids(refpath, upper_leaflet_ref_query, residue, ref_atoms[residue])

    if leaflet == "upper":
        leaflet_rsqs = upperleaflet_rsq
    elif leaflet == "lower":
        leaflet_rsqs = lowerleaflet_rsq
    else:
        print(f"error: {leaflet} not a leaflet")

    #get list of mdtraj residues of the molecule of interest
    rlist = [r for r in trj.top.residues if (r.name == residue and int(str(r)[4:]) in leaflet_rsqs)]

    #get lists of head atom indices, resnames, and resseqs 
    #that are all derived from rlist and thus all in the same order
    head_indices = [[[a.index for a in ri.atoms][charge_center_ind[residue]]] for ri in rlist]
    resns = [str(r)[0:4] for r in rlist] #the entries in this are all the same so it's sort of unnecessary but is useful in some downstream use cases
    resseqs = [str(r)[4:] for r in rlist]

    return head_indices, resns, resseqs


def closest_approaches(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/closest_approaches"
    stride = 1

    #load whole trajectory once
    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")
    
    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices=ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    dist_threshold = 1.3
    #trj_features = []

    #set dummy atom coordinates to protein center of mass
    all_prot_indices = trj.top.select(f"protein and name CA")
    prot_com = md.compute_center_of_mass(trj.atom_slice(all_prot_indices))
    for t in range(trj.n_frames):
        for k in range(3):
            trj.xyz[t][-1][k] = prot_com[t][k]
            

    for residue in ["POPC", "ARAN"]:
        for leaflet in ["upper", "lower"]:

            dist_trjs_all = []

            head_indices, resns, resseqs = get_residue_info(trj, protein, top_path, residue, leaflet)
            
            rsq_path = f"{upperpath}/{savefilename}-{residue}-{leaflet}-resseqs.npy"
            if os.path.exists(rsq_path):
                print(f"skipping {rsq_path}")
                continue
            #save out resseqs in an order matching the distance outputs
            np.save(rsq_path, resseqs)
            
            #featurize trajectory for each arachidonic acid and POPC
            for head_inds, resseq, resn in zip(head_indices, resseqs, resns):

                #get distance to protein center of mass over time
                dist_trj = get_distance_matrix(trj, head_inds, [trj.top.n_atoms-1])
                dist_trjs_all.append(dist_trj)

                t=np.argmin(dist_trj)
                if dist_trj[t] < dist_threshold:
                    
                    savefn = f"{upperpath}/{savefilename}-{resn}-rseq{resseq}-{leaflet}-frame{t*stride}-{str(round(dist_trj[t], 2))}nm.pdb"

                    trj[t].save_pdb(savefn)

                    print(f"frame = {t}; resseq = {resseq}")


            np.save(f"{upperpath}/{savefilename}-{residue}-{leaflet}.npy", np.stack(dist_trjs_all))

#test comment

compute_observable_on_trjs_x01_v2.compute_observable_x01(closest_approaches)
