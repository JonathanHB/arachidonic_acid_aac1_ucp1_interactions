import numpy as np
import mdtraj as md
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
import importlib
import os
import sys  

sys.path.insert(1, f'{os.getcwd()}/../utility')
import count_lipids
import compute_observable_on_trjs_x01_v3


##########################################################################################################################
#math to correct for box corners and volume fluctuations, allowing the RDF to be calculated out to a greater distance

def area(r1, r2, r_apothem):
    r_avg = (r1+r2)/2
    if r_avg <= r_apothem:
        return np.pi*(r2**2-r1**2)
    elif r_avg < r_apothem*np.sqrt(2):
        theta = np.arccos(r_apothem/r_avg)
        return np.pi*(r2**2-r1**2)*(1-4*theta/np.pi)
    else:
        return 0
    #TODO: integrate perimeters from r1 to r2 for more accurate areas at the corners and as the radius reaches the apothem
    
def mean_bin_area(xy_len, r2, r1):
    #box length
    box_p, box_x = np.histogram(xy_len, bins = 100, range = (8.5, 10), density=True)
    bincenters_apothem = (box_x[1:]+box_x[:-1])/4

    area_weighted = 0
    for apothem, p in zip(bincenters_apothem, box_p):
        area_weighted += area(r1, r2, apothem)*p

    return area_weighted/np.sum(box_p)


def calc_rdf(box_xy_len, molecule_radii, box_apothem, nlipids, plot = False, r_hc=[], lipid_name=""):

    molecule_radii = molecule_radii.flatten()

    #box_apothem = np.mean(box_xy_len)/2
    box_diagonal = box_apothem*np.sqrt(2)

    #calculate density
    molecules_per_bin = np.histogram(molecule_radii, bins=(200), range = (0,box_diagonal), density=False)
    bincenters = (molecules_per_bin[1][1:]+molecules_per_bin[1][:-1])/2

    rdf = []

    for i, hd0 in enumerate(molecules_per_bin[0]):
        
        #calculate the mean bin xy area within the periodic box centered on the protein
        # given the bin radius and box size distribution
        r1 = molecules_per_bin[1][i]
        r2 = molecules_per_bin[1][i+1]
        a = mean_bin_area(box_xy_len, r2, r1)

        #calculate the density of lipids per unit area in the bin
        rdf.append(nlipids*hd0/(a*len(molecule_radii)))

    #plot rdf
    if plot:
        plt.plot(bincenters, rdf)

        for r in r_hc:
            plt.axvline(r, color='red', linestyle='--')

        plt.axvline(box_apothem, color='grey', linestyle='--')

        plt.xlim(0,box_diagonal)
        #plt.ylim(0,1.6)

        plt.xlabel("radius (nm)")
        plt.ylabel(f"{lipid_name} per square nm")

    return bincenters, rdf


##########################################################################################################################
#query generation

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


#get coordinates of C-terminal ends of helices
#(last 4 residues of each helix)
def coords_c_end(protein, trj, dummy_ind):
    
    resseq_c_end = {
        "aac1": [9,94-4,114,197-4,211,294-4],
        "ucp1": [14,98-4,114,197-4,213,291-4]}
    
    mean_r = []
    mean_xy = []
    for r in resseq_c_end[protein]:
        resquery = " or ".join([f"resSeq {i}" for i in range(r, r+4)])
        atom_inds = trj.top.select(f"protein and name CA and ({resquery})")
        atom_coords = trj.xyz[0, atom_inds, 0:2]
        mean_xy.append(np.mean(atom_coords, axis=0))

        index_pairs = list(itertools.product(dummy_ind, atom_inds))
        mean_r.append(np.mean(md.compute_distances(trj, index_pairs, periodic=True, opt=True).reshape((trj.n_frames, len(atom_inds))), axis = 1))

    return np.stack(mean_xy), np.stack(mean_r)


def get_query(protein, residue, leaflet, refpath):
    #parameters
    names = {"ARAN":"C1", "POPC":"P"}
    name = names[residue]

    #select lipids or fatty acids
    if residue == "POPC":
        init_rsqs = {"aac1":1, "ucp1":9}
        upper_leaflet_ref_query = f"resSeq {init_rsqs[protein]} and name CA"
        upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet = count_lipids.count_lipids(refpath, upper_leaflet_ref_query, residue, name)
        if leaflet == "upper":
            leaflet_rsq = upperleaflet_rsq
            nleaflet = upperleaflet
        elif leaflet == "lower":
            leaflet_rsq = lowerleaflet_rsq
            nleaflet = lowerleaflet
        else: 
            raise ValueError("leaflet must be upper or lower")
        
        c_l_query = " or ".join(["resSeq "+str(i) for i in leaflet_rsq])
        ligand_query = f"resname {residue} and name {name} and ({c_l_query})"

    elif residue == "ARAN":
        c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
        invstring = ""
        if leaflet == "lower":
            invstring = " not "
        c_aran_query = " or ".join(["resSeq "+str(i) for i in c_aran[protein]])
        ligand_query = f"resname {residue} and name {name} and {invstring}({c_aran_query})"
        nleaflet = 8 #number of ARAN in each leaflet

    else:
        raise ValueError("residue must be POPC or ARAN")

    return name, f"({ligand_query})", nleaflet


##########################################################################################################################
#calculate RDFs

def calculate_radial_distances(trj, lipid_query, protein, plot = False):

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
    lipid_atoms = trj_xy.top.select(lipid_query)

    #calculate pairwise distances
    # doing it this way (i.e. using MDtraj's minimum image distance calculation for the phosphates) 
    # is important because reimaging centers molecules so the phosphate may be over the periodic boundary
    index_pairs = list(itertools.product(dummy_ind, lipid_atoms))
    lipid_radii = md.compute_distances(trj_xy, index_pairs, periodic=True, opt=True).reshape((trj_xy.n_frames, len(lipid_atoms))) #probably needlessly slow


    #####################################################################################################################
    #functions below are peripheral to this function's core purpose but are included here because it is convenient
    #####################################################################################################################

    #---------------------------part 3: calculate radii of C ends of transmembrane helices-------------------------------
    #this is included in this function because the arguments are already conveniently defined 
    helix_xy, helix_radii = coords_c_end(protein, trj_xy, dummy_ind)
    
    #---------------------------part 4: distribution of box x lengths-------------------------------
    #this is the same as the y lengths with a semiisotropic barostat
    #this is only done in this function because it is used if plot=True
    box_xy = trj.unitcell_lengths[:,0]

    #---------------------------part 5: visualize lipid distribution as a sanity check (returns nothing)-------------------------------
    #this is included in this function because the arguments are already conveniently defined 
    # and there's not any other obvious context in which you would use it
    # if plot:
    #     plot_lipid_xy(trj_xy, lipid_atoms, box_xy, prot_com)

    return lipid_radii, helix_radii, box_xy


def processing_part_1(trj_path, top_path, protein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/rdfs/"
    _savefilename = upperpath + savefilename

    residues = ["ARAN", "POPC"]
    leaflets = ["upper", "lower"]

    dummy_atom_query = {"aac1":"resSeq 590", "ucp1":"resSeq 576"}

    #loop over proteins, residues, and leaflets
    queries = []
    for residue in residues:
        for leaflet in leaflets:

            name, ligand_query, nlipids = get_query(protein, residue, leaflet, top_path)
            queries.append(ligand_query)

    #load trajectories
    ref_query = f"(protein and name CA) or {' or '.join(queries)} or {dummy_atom_query[protein]}"

    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select(ref_query)
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    for residue in residues:
        for leaflet in leaflets:

            resname, ligand_query, nlipids = get_query(protein, residue, leaflet, top_path)        
            r_all, r_hc, box_xy = calculate_radial_distances(trj, ligand_query, protein, plot=False)
            np.save(_savefilename+f"-{residue}-{leaflet}-radii.npy", r_all)
            np.save(_savefilename+f"-{residue}-{leaflet}-hc.npy", r_hc)
            np.save(_savefilename+f"-{residue}-{leaflet}-boxxy.npy", box_xy)

#UNCOMMENT FOR THE FIRST PART OF PROCESSING
#compute_observable_on_trjs_x01_v3.compute_observable_x01(processing_part_1)


def processing_part_2():

    toppath_upper = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories"
    inputpath     = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories_links"
    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/rdfs/"

    servers = ["wynton", "degrabo"]
    proteins = ["aac1", "ucp1"]
    residues = ["ARAN", "POPC"]
    leaflets = ["upper", "lower"]

    #loop over servers, proteins, and parallel runs
    for server in servers:
        for protein in proteins:
            toppath = f"{toppath_upper}/{server}/{protein}/input/seg_0035.gro"
            for run in range(1,5):
                #not all servers have the same number of parallel simulations
                if run > 2 and server == "degrabo": continue

                for residue in residues:
                    for leaflet in leaflets:

                        name, ligand_query, nlipids = get_query(protein, residue, leaflet, toppath)

                        radii = []
                        helix_radii = []
                        box_xy = []

                        for seg in range(1,20):
                            segstr = str(seg).zfill(2)

                            trjpath = f"{inputpath}/{server}/{protein}/run0{run}/seg{segstr}.xtc"
                            if not os.path.exists(trjpath): break
                            
                            #upper level directory information is supplied by the calc_obs function
                            _savefilename = upperpath+f"{server}/{protein}/run0{run}/{protein}-{server}-run0{run}-seg{segstr}"
                            
                            radii.append(np.load(_savefilename+f"-{residue}-{leaflet}-radii.npy"))
                            helix_radii.append(np.load(_savefilename+f"-{residue}-{leaflet}-hc.npy"))
                            box_xy.append(np.load(_savefilename+f"-{residue}-{leaflet}-boxxy.npy"))

                        #print(radii[0].shape)
                        #print(helix_radii)
                        #print(box_xy[0].shape)
                        r = np.concatenate(radii)
                        r_hc = np.concatenate(helix_radii, axis=1).transpose()
                        box_xy = np.concatenate(box_xy)

                        #print(r_hc.shape)

                        print(f"{len(r)} frames")

                        box_apothem = 4.5 #np.mean(box_xy)/2 #this ensures that all RDFs are calculated with the same bins

                        bincenters_all = []
                        rdf_all = []
                        hc_avg_all = []

                        increment = 5000
                        for fmax in range(increment, len(box_xy), increment):

                            bincenters, rdf = calc_rdf(box_xy[:fmax], r[:fmax], box_apothem, nlipids, plot = True, r_hc=[], lipid_name="")
                            bincenters_all.append(bincenters)
                            rdf_all.append(rdf)
                            
                            hc_avg = np.mean(r_hc[:fmax], axis = 0)
                            hc_avg_all.append(hc_avg)

                        _savefilename = upperpath+f"outputs/{server}/{protein}/run0{run}/{protein}-{server}-run0{run}"
                        np.save(_savefilename+f"-{residue}-{leaflet}-bincenters.npy", np.stack(bincenters_all))
                        np.save(_savefilename+f"-{residue}-{leaflet}-rdf.npy", np.stack(rdf_all))
                        np.save(_savefilename+f"-{residue}-{leaflet}-helix-radii.npy", np.stack(hc_avg_all))


processing_part_2()
