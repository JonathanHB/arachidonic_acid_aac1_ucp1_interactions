# CELL 1
# IMPORTS
##########################################################################################

#general
import time
import os
import sys

import itertools
import numpy as np

import matplotlib.pyplot as plt

#md-specific
import mdtraj as md

#my functions
sys.path.insert(1, f'{os.getcwd()}/../utility')

# import load_all_trjs
# import load_anatale_trjs
import count_lipids


# CELL 2
# RESIDUE INFORMATION FOR BINDING SITE CLASSIFICATION
##########################################################################################

# def helix_resqueries(protein):

#     terminal_residues = {
#         "aac1":
#        [[9,   36],
#         [74,  94],
#         [114, 141],
#         [177, 197],
#         [211, 238],
#         [274, 294]],
#         "ucp1":
#        [[14,  41],
#         [78,  98],
#         [114, 141],
#         [177, 197],
#         [213, 240],
#         [271, 291]]
#     }
    
#     resqueries = []
#     for rt in terminal_residues[protein]:
#         helix_resseqs = [i for i in range(rt[0], rt[1]+1)]
#         #print("color red, resi " + "+".join([str(i) for i in helix_resseqs]))
#         resqueries.append(" or ".join([f"resSeq {i}" for i in helix_resseqs]))
    
#     return resqueries


def site_def_resseqs(protein):

    terminal_residues = {
        "aac1":
       [[5,   27],
        [94,  74],
        [109, 132],
        [200, 177],
        [209, 229],
        [294, 274]],
        #deeper symmetrical version
    #    [[9,   27],
    #     [94,  74],
    #     [114, 132],
    #     [197, 177],
    #     [211, 229],
    #     [294, 274]],
        "ucp1":
       [[14,  30],
        [78,  98],
        [114, 130],
        [177, 197],
        [213, 229],
        [271, 291]]
    }

    return terminal_residues[protein]


# CELL 3
# GENERATE MDTRAJ QUERY FOR SPECIFIED (RESIDUE, LEAFLET, PROTEIN) COMBINATION
##########################################################################################

def get_ligand_query(protein, residue, leaflet, refpath):
    #parameters
    names = {"ARAN":["C1", "C5", "C10", "C15", "C20"], "POPC":["P"]}
    namequery = " or ".join(["name "+ name for name in names[residue]])

    #select lipids or fatty acids
    if residue == "POPC":
        init_rsqs = {"aac1":1, "ucp1":9}
        upper_leaflet_ref_query = f"resSeq {init_rsqs[protein]} and name CA"
        upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet = count_lipids.count_lipids(refpath, upper_leaflet_ref_query, residue, names[residue][0])
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



# CELL 6
# MATHEMATICS OF FACE INTERSECTION CALCULATION USING TRIPLE PRODUCTS
##########################################################################################
#see https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d


def sgnv(a,b,c,d):
    return np.sum(np.multiply(d-a, np.cross(b-a, c-a)), axis = 1)
    
    #triple product
    #this does not work with array inputs but illustrates what happens to each layer of the input arrays
    #return np.dot(d-a, np.cross(b-a, c-a)) 


#takes 2d array inputs
def facecrossings(p1,p2,p3, q1,q2):
    
    #check if q1 and q2 are on opposite sides of triangle p1p2p3
    opp_sgnvs = np.multiply(sgnv(p1, p2, p3, q1), sgnv(p1, p2, p3, q2))
    opp_sides = np.where(opp_sgnvs<0, 1, 0)

    #check if the line containing q1 and q2 intersects triangle p1p2p3
    lit_sgnvs = np.stack([sgnv(q1, q2, p1, p2), sgnv(q1, q2, p2, p3), sgnv(q1, q2, p3, p1)])
    #print(lit_sgnvs.shape)

    plus_signs = np.where(lit_sgnvs>0, 1, 0)
    minus_signs = np.where(lit_sgnvs<0, 1, 0)
    all_same_sign = np.product(plus_signs, axis=0) + np.product(minus_signs, axis=0)

    #return 1 only for frames where q1q2 passes through p1p2p3
    return np.multiply(opp_sides, all_same_sign)


# CELL 7
# FEATURIZE TRAJECTORY TO IDENTIFY BOUND AND UNBOUND STATES AND BINDING SITES
##########################################################################################

#returns a (7 x n_frames) array of 1s and 0s describing whether there is a AA bound in groove i on frame j. 
# the last index is for FA entirely within the cavity

def binding_state(trj, dummy_atom_query, protein):
    
    #------------------------------------------------------------------------------------
    # PROTEIN CENTER OF MASS

    #set dummy atom coordinates to protein center of mass
    all_prot_indices = trj.top.select(f"protein and name CA")
    dummy_ind = trj.top.select(dummy_atom_query)[0]
    prot_com = md.compute_center_of_mass(trj.atom_slice(all_prot_indices))
    for t in range(trj.n_frames):
        for k in range(3):
            trj.xyz[t][dummy_ind][k] = prot_com[t][k]


    #------------------------------------------------------------------------------------
    # DISTANCES FROM PROTEIN CENTER
    #calculate with atoms are within threshold of the protein COM (a n_frames x n_AA array of 1s and 0s)

    #protein and FA atom selections
    atoms = ["C1", "C5", "C10", "C15", "C20"]
    c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
    helix_ends = site_def_resseqs(protein)
    helix_ends_inds = []
    for h in helix_ends:
        helix_ends_inds.append([trj.top.select(f"resSeq {h[0]}")[0], trj.top.select(f"resSeq {h[1]}")[0]])

    #the code below will crash with an uninformative error message if the wrong number of brackets are used
    events_by_helix_fa = [[],[],[],[],[],[],[]]
    events_by_helix = [[],[],[],[],[],[],[]]

    #save examples of FA bound in each site
    save_examples = False
    helices_written = []

    plot_binding_trjs = False

    #keeps track of whether FA could be in the cavity given which side of the plane between each pair of helices they're on
    cavity_bound_all = []

    #loop over pairs of consective helices
    for hi in range(6): 
        print(f"groove {hi+1}")

        cavity_bound_r = [] #cavity binding information for all FA and the current helix

        #loop over FA molecules
        for r in c_aran[protein]:

            #an array counting how many pairs of FA atoms cross the planes defining the binding site in each frame
            # it is later digitized since usually only one or 0 carbons cross the plane
            ingroove = np.zeros(trj.n_frames)

            #loop over pairs of consecutive FA atoms
            for ai in range(len(atoms)-1): 

                #select two atoms
                a1j = trj.top.select(f"resname ARAN and resSeq {r} and name {atoms[ai]}")[0]
                a2j = trj.top.select(f"resname ARAN and resSeq {r} and name {atoms[ai+1]}")[0]

                #calculate whether the line between the atoms intersects each binding site plane
                #  
                triangle_1_intersection = facecrossings(trj.xyz[:, helix_ends_inds[hi][0]], trj.xyz[:, helix_ends_inds[hi][1]],       trj.xyz[:, helix_ends_inds[(hi+1)%6][1]], trj.xyz[:, a1j], trj.xyz[:, a2j])
                triangle_2_intersection = facecrossings(trj.xyz[:, helix_ends_inds[hi][0]], trj.xyz[:, helix_ends_inds[(hi+1)%6][0]], trj.xyz[:, helix_ends_inds[(hi+1)%6][1]], trj.xyz[:, a1j], trj.xyz[:, a2j])

                ingroove += triangle_1_intersection + triangle_2_intersection

                # #example saving not tested recently, may contain bugs
                # if save_examples:
                #     if (hi not in helices_written): #write one example per site
                #         for t, tr in enumerate(triangle_1_intersection + triangle_2_intersection):
                #             if tr == 1: #look for a frame in which the current pair of atoms is bound. If none exists nothing happens.
                #                 pdb_path = f"/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/visualization/{protein}"
                #                 c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
                #                 trj[t].save_pdb(f"{pdb_path}/seg_{mtd[0][-12:-6]}_{mtd[1]}_AA_site_{hi+1}_resi_{r}_atom_{atoms[ai]}_{t*10}_v2.pdb")
                #                 print(f"saved site {hi+1} resi {r} atom {atoms[ai]} frame {t*10}")
                                
                #                 helices_written.append(hi)
                #                 break
            
            #record whether any pair of atoms spans the binding site plane for the current FA and site
            ingroove_bin = np.where(ingroove>0, 1, 0)
            events_by_helix_fa[hi].append(ingroove_bin)

            #plot the binding trajectory of current FA
            if plot_binding_trjs:
                plt.plot(ingroove_bin)
                plt.legend([str(i) for i in c_aran[protein]])
            
            #check if all atoms of the current FA are on the inside side of the current site's reference planes
            sgnv_interior = []
            for a in atoms:
                ai = trj.top.select(f"resname ARAN and resSeq {r} and name {a}")[0]

                triangle_1_sgnv = sgnv(trj.xyz[:, helix_ends_inds[hi][0]], trj.xyz[:, helix_ends_inds[hi][1]],       trj.xyz[:, helix_ends_inds[(hi+1)%6][1]], trj.xyz[:, ai])
                t1_sv_bin = np.where(triangle_1_sgnv<0, 1, 0)
                sgnv_interior.append(t1_sv_bin)

                triangle_2_sgnv = sgnv(trj.xyz[:, helix_ends_inds[hi][0]], trj.xyz[:, helix_ends_inds[(hi+1)%6][1]], trj.xyz[:, helix_ends_inds[(hi+1)%6][0]], trj.xyz[:, ai])                
                t2_sv_bin = np.where(triangle_2_sgnv<0, 1, 0)
                sgnv_interior.append(t2_sv_bin)

            cavity_bound = np.where(np.sum(np.stack(sgnv_interior), axis = 0) == 10, 1, 0)
            cavity_bound_r.append(cavity_bound)

        #check if any FA is bound in the current groove at each time point
        events_by_helix_fa[hi] = np.stack(events_by_helix_fa[hi])
        events_by_helix[hi] = np.where(np.max(np.stack(events_by_helix_fa[hi]), axis = 0) == 1, 1, 0)

        #add potential cavity binding data for current helix
        cavity_bound_all.append(np.stack(cavity_bound_r))

        #show plot of each FA's binding in the current site        
        plt.show()

    #find FA which meet the cavity binding criteria for the binding planes for all 6 sites
    cavity_trj_byfa = np.where(np.sum(np.stack(cavity_bound_all), axis = 0) == 6, 1, 0)

    #plot which FA are bound in the C-side cavity as a function of time
    if plot_binding_trjs:
        print("cavity binding")
        for i in range(len(cavity_trj_byfa)):
            plt.plot(cavity_trj_byfa[i])

        plt.legend([str(i) for i in c_aran[protein]])
        plt.show()

    #find frames in which any FA is bound in the central cavity and add it to the list specifying which interhelical sites have bound FA
    events_by_helix_fa[6] = np.stack(cavity_trj_byfa)
    events_by_helix[6] = np.where(np.max(np.stack(cavity_trj_byfa), axis = 0) == 1, 1, 0)

    #convert to array and return
    events_by_helix = np.stack(events_by_helix)
    events_by_helix_fa = np.stack(events_by_helix_fa)

    return events_by_helix, events_by_helix_fa
