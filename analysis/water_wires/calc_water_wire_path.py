import numpy as np
import itertools

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse.csgraph import shortest_path
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt

import mdtraj as md

import time
import os
import sys

sys.path.insert(1, f'{os.getcwd()}/../utility')
import compute_observable_on_trjs_x01_v3 

#Jonathan Borowsky
#Grabe lab
#8/27/25

#for centered trajectories periodicity checks are not needed
#THIS CODE IS ONLY RELIABLE FOR CENTERED TRAJECTORIES. PLEASE CENTER YOUR TRAJECTORY USING GROMACS.
def get_water_cylinder_inds(frame_xyz_prot, water_o):
    #get the indices of all water oxygens in a cylinder around the protein
    #note that this does not work across PBC

    #extent of cylinder from protein COM in nm
    z_dist_p = 3.5 #6
    z_dist_n = 4 #-2.5
    r_dist = 2
    r_sqdist = r_dist**2
    
    #restrict water which can be used for path endpoints 
    # to those at least [pad] from the edge of the cylinder. 
    # As an endpoint on the surface of the cylinder can only be approached from one side, 
    # the average best path to it is worse than the average best path through a water within the cylinder
    # this effect shows up noticeably on plots of water paths and affects the value of the progress coordinate 
    # once the protein channel is fully or mostly open
    # we therefore calculate a smaller cylinder and require the endpoints to be inside of it
    pad = 0.4
    
    z_dist_p_end = z_dist_p - pad
    z_dist_n_end = z_dist_n - pad
    r_dist_end = r_dist - pad
    r_sqdist_end = r_dist_end**2


    #squared radii of the water oxygens
    sqradii = np.sum(frame_xyz_prot[water_o][:, :2]**2, axis=1)

    #find water oxygens in outer cylinder
    in_radius = np.where(sqradii < r_sqdist, 1, 0)
    in_z = np.where((frame_xyz_prot[water_o][:, 2] < z_dist_p) & (frame_xyz_prot[water_o][:, 2] > -z_dist_n), 1, 0)
    cyl_o_inds_inds = np.nonzero(np.multiply(in_radius, in_z)==1)
    cyl_o_inds = list(water_o[cyl_o_inds_inds])

    #find water oxygens in inner cylinder for endpoints
    in_radius_end = np.where(sqradii < r_sqdist_end, 1, 0)
    in_z_end = np.where((frame_xyz_prot[water_o][:, 2] < z_dist_p_end) & (frame_xyz_prot[water_o][:, 2] > -z_dist_n_end), 1, 0)
    cyl_o_inds_inds_end = np.nonzero(np.multiply(in_radius_end, in_z_end)==1)
    cyl_o_inds_end = list(water_o[cyl_o_inds_inds_end])

    #print(f"cylinder contains {len(cyl_o_inds)} waters")
    
    #get waters at either end of the smaller cylinder for path calculation
    cyl_o_z = list(frame_xyz_prot[cyl_o_inds_end,2])

    #get indices of the terminal waters in the cyl_o_inds array rather than the cyl_o_inds_end array
    c_water_ind = cyl_o_inds.index(cyl_o_inds_end[np.argmin(cyl_o_z)])
    m_water_ind = cyl_o_inds.index(cyl_o_inds_end[np.argmax(cyl_o_z)])

    return cyl_o_inds, c_water_ind, m_water_ind


#lack of same-residue checking is not needed for water-only wire with terminal ARAN
#for centered trajectories periodicity checks are not needed
def get_sqdist_matrix(coords_all, indices):

    coords = coords_all[indices]

    distances = squareform(pdist(coords, 'sqeuclidean')) #~1/4 runtime

    #TODO: write an option to identify connections between atoms on the same residue to exclude them 
    #since protons can't go directly from one glutamate O to the other

    #identify long connections and exclude them
    distances_zeroed = np.where(distances > 1, 0, distances)

    #convert to sparse matrix (this is actually the slowest part, accounting for over 1/4 of the program's entire runtime)
    distances_csr = csr_matrix(distances_zeroed)

    return distances, distances_csr


############################################################################################################
#                                    WATER WIRE PATH CLASSIFICATION
############################################################################################################


def get_vertex_inds(trj, protein):

    pro_resseqs ={"ucp1": [32, 78, 132, 178, 231, 272], "aac1": [27, 74, 132, 178, 229, 275]} #for we37 ucp1 numbering
    val_thr_resseqs = {"ucp1":[38, 138, 237], "aac1": [33, 138, 235]} #for we37 ucp1 numbering

    #for checking that the indices are right
    #print("+".join([str(i) for i in pro_resseqs[protein]+val_thr_resseqs[protein]]))

    val_thr_str = " or ".join([f"resSeq {r}" for r in val_thr_resseqs[protein]]) #resSeq is case sensitive
    val_thr_ca_inds = trj.top.select(f"({val_thr_str}) and name CA")

    pro_str = " or ".join([f"resSeq {r}" for r in pro_resseqs[protein]]) #resSeq is case sensitive
    pro_ca_inds = trj.top.select(f"({pro_str}) and name CA")
    
    return val_thr_ca_inds, pro_ca_inds


#determine if the line between atoms at coordinates d and e crosses a face* of the pyramid, and if so which one. 
# *one of the angles (in the geometer's sense of the word) between formed by the top of the pyramid and a pair 
# of adjacent base vertices. The base of the pyramid is not a face for these purposes.
#---------------------------------
#Parameters:
# vertices: a 7x3 array where each row is the coordinates of a vertex, starting with the top of the pyramid
# d, e: lists of 3 floats: the coordinates of 2 points
#---------------------------------
#Returns:
# an integer describing whether the line between d and e passes through the pyramid, and if so which one
#    -1: the line segment de does not intersect one of the pyramid's triangles
#    0-2 inclusive: de intersects one triangle from the corresponding pair (or both but this has never to my 
#       knowledge happened and would be geometrically ridiculous given the biological context; 
#       at least one single-triangle intersection must occur.)

#NOTE: this method actually finds if de intersects the plane containing angle bac rather than triangle bac
def facecrossings(vertices, d, e):
        
    for i in range(6):
        
        a = vertices[0]
        b = vertices[i+1]
        c = vertices[(i+1)%6 + 1]
        
        #points with respect to the point of the pyramid
        ba = b-a
        ca = c-a
        da = d-a
        ea = e-a
        
        normvec = np.cross(ba, ca)
        on_opp_side = np.dot(normvec,da)*np.dot(normvec,ea)
    
        #if one or both points lie exactly in the plane
        if on_opp_side == 0:
            #statistically this should barely ever happen
            print("warning: 0")
        
        #if the points lie on opposite sides of the plane (and the proton path thus crosses the plane)
        if on_opp_side <= 0:
           
            #get the vectors between (d. and e.) and a. in the plane abc
            da_perp = np.dot(normvec, da)/np.dot(normvec,normvec)*normvec
            da_parr = da - da_perp
            ea_perp = np.dot(normvec, ea)/np.dot(normvec,normvec)*normvec
            ea_parr = ea - ea_perp

            #the point where line de intersects triangle abc
            intersection_point = (da_parr*np.linalg.norm(ea_perp) + ea_parr*np.linalg.norm(da_perp))/(np.linalg.norm(da_perp) + np.linalg.norm(ea_perp))

            #negative if b and c are on opposite sides of the line between a and the intersection
            sgnbc = np.dot(np.cross(intersection_point, ba), np.cross(intersection_point, ca))
            
            #negative if a and c are on opposite sides of the line between b and the intersection
            #sgnca = np.dot(np.cross(intersection_point-ba, -ba), np.cross(intersection_point-ba, c-b))
        
            #angles (wrt a) between the intersection and points b and c
            i_b_angle = np.arccos(np.dot(intersection_point, ba)/(np.linalg.norm(intersection_point)*np.linalg.norm(ba)))
            i_c_angle = np.arccos(np.dot(intersection_point, ca)/(np.linalg.norm(intersection_point)*np.linalg.norm(ca)))
            total_angle = i_b_angle+i_c_angle
            
            if sgnbc<=0 and total_angle < np.pi: #and sgnca<=0: #see notebook
                if i == 0 or i == 1:
                    return 0
                elif i == 2 or i == 3:
                    return 1
                else:
                    return 2
    
    return -1


# take a water pathway and use facecrossings to figure out which third of the protein it passes through
#---------------------------------
#Parameters:
# trj: an MDtraj trajectory (only the last frame is used)
# path_o_inds: list of strings: indices of water wire forming atoms in mdtraj frame
# vertices: a 7x3 array where each row is the coordinates of a vertex, starting with the top of the pyramid
#---------------------------------
#Returns:
# an integer describing which third of the protein the water wire passes through, or an edge case
#    0-2 inclusive: the corresponding third of the protein
#    3: the water wire was never observed to cross the planes defined by vertices. 
#        This indicates a bug or badly selected water wire endpoints
#    4: the water wire made multiple crossings through different thirds of the protein

def classify_path(xyz, path_o_inds, vertices):

    
    #check each pair of points along the water wire to find the one crossing the pyramid
    fcs = []
    for i in range(len(path_o_inds)-1):
        fc = facecrossings(vertices, xyz[path_o_inds[i]], xyz[path_o_inds[i+1]])
        if fc != -1:
            fcs.append(fc)
    
    #return the face that the water wire passed through
    if len(fcs) == 1:
        return fcs[0]
        
    #if the water wire looped back and crossed multiple times
    elif len(fcs)>1:
        #If all crossings passed through the same third of the protein, 
        #report this as if it were a regular water wire through that part of the protein
        #from what I've seen this occurs with multiple-water-wide water wires 
        #through which the minimum spanning tree takes circuitous paths
        #but these are not fundamentally different from slightly less convoluted wires 
        #through the same channels and do not seem to involve independent return loops
        fc0 = fcs[0]
        for fc in fcs:
            if fc != fc0:
                #if the water wire passes through multiple parts of the protein, return a nonstandard PC value (4)
                return 4
        return fc0
        
    #report cases where the path cannot be found or classified
    else:
        #note that this print statement can and should crash westpa
        #it has fortunately never been observed (except as a failure mode in the case of erroneous inputs)
        print("path not classified and/or did not pass through any face")
        return 3


############################################################################################################

#TODO can we remove frame as an argument to make this mdtraj independent?
def get_water_pc(centered_coords, water_o, vertices, printcommands = False):

    #t1 = time.time()

    #get cylinder waters
    cyl_o_inds, c_water_ind, m_water_ind = get_water_cylinder_inds(centered_coords, water_o)

    if printcommands:
        #indices printed here do not play nicely with pdb files for some reason
        #print([i%3 for i in cyl_o_inds])
        print(f"hide spheres; show spheres, index {'+'.join([str(i+1) for i in cyl_o_inds])}")

    #t2 = time.time()
    #print(f"cylinder waters identified in {t2-t1} seconds")

    #calculate distances between the waters in the cylinder
    o_dist_mat, o_dist_mat_csr = get_sqdist_matrix(centered_coords, cyl_o_inds)

    #t3 = time.time()
    #print(f"distances calculated in {t3-t2}")

    #calculate the minimum spanning tree, which contains the minmax path between every pair of waters
    #mst is some kind of scipy sparse matrix or graph object
    mst = minimum_spanning_tree(o_dist_mat_csr)

    #distance matrix containing the mst, with 0 used for nonexistent edges
    #note that this matrix is not symmetrical; each nonzero entry has a 0 entry at the transposed position
    #mst_arr = mst.toarray()

    #t4 = time.time()
    #print(f"minimum spanning tree calculated in {t4-t3} seconds")

    #determine the shortest paths between the water with the smallest z coordinate and the water with the largest 
    dist_matrix, predecessors = shortest_path(mst, directed=False, indices=c_water_ind, return_predecessors=True)

    #t5 = time.time()
    #print(f"shortest path calculated in {t5-t4} seconds")

    #traverse the predecessor list to find the shortest path
    nodelist = [m_water_ind]
    path_dists = []

    for t in range(len(predecessors)):
        #append distance from the previous to current node
        if t != 0:
            path_dists.append(o_dist_mat[nodelist[-2]][nodelist[-1]])

        if predecessors[nodelist[-1]] < 0:
            break
        #add the next node to the path
        nodelist.append(predecessors[nodelist[-1]])

    nodelist_structure_inds = np.array(cyl_o_inds)[np.array(nodelist)]

    #water wire coordinates
    xyz = centered_coords[nodelist_structure_inds]
    
    if printcommands: 
        print(f"hide spheres; show spheres, index {'+'.join([str(n+1) for n in nodelist_structure_inds])}")

    #t6 = time.time()
    #print(f"minmax MST edge length took {t6-t5} seconds")

    #print("Total time: ", t7-t1)

    #if the minmax water wire gap is longer than 1 nm, 
    # the water adjacency graph is not connected due to the sparse matrix implementation above, 
    # and the path distances are nonsense
    if len(path_dists) != 0:
        #only calculate the path class if there is a connected MST such that the path can be classified  
        path_class = classify_path(centered_coords, nodelist_structure_inds, vertices)
        return path_dists, xyz, nodelist_structure_inds, path_class
    else:
        #5 is a dummy path class for cases where the gap is so large that the MST is disconnected
        return [1], xyz, nodelist_structure_inds, 5



#---------------------------------------------------------------------------------------
#MAIN METHOD
    
def save_wire_pc(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/water_wires_no_aa/"
    _savefilename = upperpath + savefilename

    queries = ["resname HOH and element O"]#,
               #"resname ARAN and (name O1 or name O2)"]
    threshold = 0.33 #nm; distance below which water wires are considered connected

    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    #trj.save_gro("/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/water_wires/degrabo/aac1/run02/aac1-degrabo-run02-segnum15-1010-0.29nm.gro")

    print(f"trajectory has {trj.n_frames} frames")

    #assemble queries
    combined_queries = " or ".join([f"({q})" for q in queries])

    check_queries = True
    if check_queries:
        for q in queries:
            if len(trj.top.select(q)) == 0:
                print(f"error: no results for {q}")

    water_o = trj.top.select(combined_queries)

    val_thr_ca_inds, pro_ca_inds = get_vertex_inds(trj, protein)

    #coordinates relative to protein center of mass, not PBC adjusted
    prot_com = md.compute_center_of_mass(trj, select = "protein")
    trj_xyz_prot = trj.xyz - np.tile(prot_com, (trj.xyz.shape[1],1,1)).transpose(1,0,2)


    pcs = []
    paths = []

    for fi in range(0, trj.n_frames):
        if fi%50 == 0:
            print(fi)

        #t1 = time.time()

        #average valine and threonine CAs to get 'top' of pyramid
        top_vertex = np.mean(trj_xyz_prot[fi, val_thr_ca_inds], axis=0)
        base_vertices = trj_xyz_prot[fi, pro_ca_inds]
        vertices = np.vstack([top_vertex, base_vertices])

        #t11 = time.time()
        path_dists, xyz, nodelist, path_class = get_water_pc(trj_xyz_prot[fi], water_o, vertices, printcommands = False)
        #print(f"main time = {t22-t11}")
        # print(f"frame {fi} took {t2-t1} seconds")

        paths.append(path_class)

        maxdist = np.sqrt(max(path_dists))
        pcs.append(maxdist)

        #t2 = time.time()
        #print(f"loop time = {t2-t1}")

        #save coordinates of water wire atoms
        if maxdist <= threshold:
            np.save(f"{_savefilename}-frameind{fi}-xyz.npy", xyz)

        if path_class == 3:
            trj[fi].save_gro(f"{_savefilename}-frameind{fi}-pathclass3")

    #save water wire lengths and paths
    np.save(_savefilename+"-gaps" + ".npy", pcs)

    #print(paths)
    paths = np.array(paths).astype("uint8")
    np.save(_savefilename+"-paths" + ".npy", paths)

    return pcs, paths



compute_observable_on_trjs_x01_v3.compute_observable_x01(save_wire_pc)


