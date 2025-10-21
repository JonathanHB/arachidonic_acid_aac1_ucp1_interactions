import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import mdtraj as md

#####################################################################################
#-------------------------------------pdb file io------------------------------------

#pdb reader
#pdb file format: https://www.wwpdb.org/documentation/file-format
def pdb_loader(input_fn, output_fn, atoms, bfactors):

    #format: [[index, atom number, atom name, x, y, z, whether the atom is hydrogen],...]

    x = 0
    with open(output_fn, "w") as fo:    
        with open(input_fn, "r") as fi:
            for line in fi:
                if line[0:4] == "ATOM":

                    atomstring = line[17:20].strip() + str(int(line[22:26].strip())) + "-" + line[12:16].strip()
                    
                    if line[17:20].strip() == "HSD":
                        atomstring = "HIS" + str(int(line[22:26].strip())) +"-"+ line[12:16].strip()
                    elif line[17:20].strip() == "ILE" and line[12:16].strip() == "CD":
                        atomstring = line[17:20].strip() + str(int(line[22:26].strip())) +"-"+ "CD1"
                                            
                    if atomstring in atoms:
                        fo.write(line[0:61] + str(round(bfactors[atoms.index(atomstring)], 2)).rjust(5) + line[66:])
                    else:
                        fo.write(line)


def save_pdb_bfactors(gro_path, pdb_path, output_path, val_by_resi, atom_query, suffix=""):

    #----------------------------------------------------------
    #get heavy atom names to go with indices
    gro_frame = md.load(gro_path)
    
    gro_indices = gro_frame.top.select(atom_query)
    fta = [a for a in gro_frame.top.atoms]
    heavy_atoms = [str(fta[i]) for i in gro_indices]
    
    absmax_r = max([abs(i) for i in val_by_resi])
    
    bfactors = [ri/absmax_r for ri in val_by_resi]

    pdb_loader(pdb_path, output_path +"/"+ suffix + ".pdb", heavy_atoms, bfactors) 

#-------------------------------------end pdb file io--------------------------------
#####################################################################################


#list of residues by helix
def helix_resqueries(protein):

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
    
    #rq_inds = [np.array([i for i in range(0,rt[1]+1-rt[0])]) for rt in terminal_residues[protein]]
    resseqs = []
    resqueries = []
    for rt in terminal_residues[protein]:
        helix_resseqs = [i for i in range(rt[0], rt[1]+1)]
        resseqs += helix_resseqs
        #print("color red, resi " + "+".join([str(i) for i in helix_resseqs]))
        resqueries.append(" or ".join([f"resSeq {i}" for i in helix_resseqs]))
    
    return resqueries, resseqs


def contact_freq_byresi(contacts_by_atom, binding_events, atoms_by_residue, gro_path, pdb_path, output_path, atom_query, protein):

    #####################################################################################################
    # code to make figure of AA contact frequency by residue for each binding site
    #####################################################################################################

    #loop over binding sites
    for si in range(7):

        #contact frequency by residue for AA bound at current site
        freqs_by_residue = []

        #loop over residues
        for abr in atoms_by_residue:
            n_contacts = 0
            #loop over atoms in residue
            for a in abr:
                n_contacts += np.sum(np.multiply(contacts_by_atom[:,:,a], binding_events[si]))

            #copilot suggested this normalization factor and I have to say it's better than what I had in mind
            freqs_by_residue.append(n_contacts / np.sum(binding_events[si])) 

        save_pdb_bfactors(gro_path, pdb_path, output_path, freqs_by_residue, atom_query, suffix=f"-{protein}-site-{si}")


#1. load protein nonhydrogen atoms and get indices of atoms in each helix
#2. load contacts
#3. 8 x n_frames x 6 binary array stating whether each AA contacts each helix at each frame
#   to do this multiply the contact array by a one hot encoding of which atoms are in the ith helix, and repeat for the other helices
#   Then look for frames where the AA crosses the plane between helices i and j and also contacts both of them

def identify_bound_states(ibs):

    #switch to x01 paths
    inputpath_main = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing" #"/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1"
    inputpath_contacts = f"{inputpath_main}/contacts"
    inputpath_planes = f"{inputpath_main}/binding_planegeom"
    toppath_upper = f"{inputpath_main}/centered_trajectories"

    servers = ["wynton", "degrabo"]
    proteins = ["aac1", "ucp1"]

    #loop over servers, proteins, and parallel runs
    for protein in proteins:

        #---------------------------get helix indices to convert atom-level contacts to helix-level contacts---------------------------
        rq, rsq = helix_resqueries(protein)

        ref_trj = md.load(f"{toppath_upper}/{servers[0]}/{protein}/input/seg_0035.gro")
        ref_inds = ref_trj.top.select("protein and not element H")
        trj = md.load(f"{toppath_upper}/{servers[0]}/{protein}/input/seg_0035.gro", atom_indices = ref_inds)

        cumulative_atoms = 0
        helix_inds = []
        for rqi in rq:
            sele = trj.top.select(rqi)
            helix_inds.append([i+cumulative_atoms for i in range(len(sele))])
            cumulative_atoms += len(sele)

        #####################################################################################################
        #convert to residue level contacts for figure added for revisions
        #####################################################################################################

        cumulative_atoms = 0
        atoms_by_residue = []
        for rseq in rsq:
            sele = trj.top.select(f"resSeq {rseq} and not element H")
            atoms_by_residue.append([i+cumulative_atoms for i in range(len(sele))])
            cumulative_atoms += len(sele)


        #-----------------------------------------------------------------------------------------------------------------------------

        contacts_all_all = []
        contacts_planes_all_all = []

        rn = 1
        for server in servers:
            for run in range(1,5):
                
                print(protein, run, server)

                planes_all = []
                contacts_all = []
                contacts_byhelix_all = []
                contacts_planes_all = []

                if server == "degrabo" and run > 2:
                    continue

                #step is in microseconds
                if server == "wynton" and (run == 3 or run == 4):
                    step = 1/10000
                else:
                    step = 1/5000

                # dists = []
 
                for seg in range(1,20):

                    file_contacts = f"{inputpath_contacts}/{server}/{protein}/run0{run}/{protein}-{server}-run{str(run).zfill(2)}-seg{str(seg).zfill(2)}-upper-ARAN-contacts-byaa.npy"
                    file_planes = f"{inputpath_planes}/{server}/{protein}/run0{run}/{protein}-{server}-run{str(run).zfill(2)}-seg{str(seg).zfill(2)}-plane-fa-crossings.npy" 
                    
                    #print(file_contacts)
                    #everything can be done in here I believe; it should probably be its own function
                    if os.path.exists(file_contacts):
                        #contacts: [n_upperleaflet_aa x n_frames x n_protein_heavy_atoms] array of 0s and 1s 
                        contacts = np.load(file_contacts).astype("uint8").transpose(2,0,1)
                        contacts_all.append(contacts)
                        #contacts.transpose(2,0,1)
                        #planes: [n_helices+1 x n_upperleaflet_aa x n_frames] 
                        planes = np.load(file_planes).astype("uint8")
                        
                        cbn, cap = ibs(contacts, planes, helix_inds)
                        #sys.exit(0)
                        planes_all.append(planes)
                        contacts_byhelix_all.append(cbn)
                        contacts_planes_all.append(cap)
                        
                # print(len(planes_all))
                planes_all = np.concatenate(planes_all, axis=2)
                contacts_all = np.concatenate(contacts_all, axis=1)
                contacts_byhelix_all = np.concatenate(contacts_byhelix_all, axis=2)
                contacts_planes_all = np.concatenate(contacts_planes_all, axis=2)
                
                contacts_all_all.append(contacts_all.astype("uint8"))
                contacts_planes_all_all.append(contacts_planes_all.astype("uint8"))

                print(planes_all.shape)
                print(contacts_all.shape)
                print(contacts_byhelix_all.shape)
                print(contacts_planes_all.shape)


                colors = ["purple", "blue", "cyan", "green", "orange", "red", "magenta", "grey"]
                legend = ["1", "2", "3", "4", "5", "6", "7", "8"]

                plt.figure(dpi=600)

                #plot interhelical binding events
                for i in range(6):
                    #plt.title(f"h{i+1}-h{(i+1)%6+1}")
                    for j in range(8):
                        #plt.plot(planes_all[i][j], color = "blue")
                        #plt.plot(contacts_all[i][j], color = "red")
                        plt.scatter([t*step for t in range(len(contacts_planes_all[i][j]))], (i+1)*contacts_planes_all[i][j]+j/32-1/8, color = colors[j], s=0.5, marker="o")

                #plot cavity binding events, which are not assigned based on contacts
                for j in range(8):
                    plt.scatter([t*step for t in range(len(planes_all[6][j]))], 7*planes_all[6][j]+j/32-1/8, color = colors[j], s=0.5, marker="o")
                
                plt.ylim(0.5,7.5)
                plt.xlim(0,7.5)
                plt.yticks([1,2,3,4,5,6,7], ["h1-h2", "h2-h3", "h3-h4", "h4-h5", "h5-h6", "h6-h1", "cavity"])
                plt.ylabel("binding site")
                plt.xlabel(r"time ($\mu$s)")
                plt.title(f"{protein} run {rn}")
                plt.legend(legend)
                plt.axvline(x = step*len(planes_all[0][0]) , color = "black", linestyle = "dashed")        

                plt.savefig(f"{inputpath_planes}/{protein}-{server}-run{str(run).zfill(2)}-binding.png")
                plt.show()

                rn += 1
                #sys.exit(0)

        #########################################################################################################
        #----------------------------------------for revisions figure--------------------------------------------

        contacts_all_all = np.concatenate(contacts_all_all, axis=1)
        contacts_planes_all_all = np.concatenate(contacts_planes_all_all, axis=2)

        helixresis_all = "or".join(rq)
        protein_query = f"protein and ({helixresis_all}) and not element H"
        
        pdb_path = md.load(f"{toppath_upper}/{servers[0]}/{protein}/input/seg_0035.pdb")
        output_path = f"{inputpath_main}/contact_freq_heatmaps"

        contact_freq_byresi(contacts_all_all, contacts_planes_all_all, atoms_by_residue, ref_trj, pdb_path, output_path, protein_query, protein)


def _identify_bound_states(contacts, planes, helix_inds):

    #convert atom-level contacts to helix level ones
    contacts_by_helix_frame_resi = []
    for hi in helix_inds:
        #print(contacts.shape)
        helix_contacts = contacts[hi]
        #print(helix_contacts.shape)
        contact_num = np.sum(helix_contacts, axis = 0)
        #ps(contact_num)
        contact_bin = np.where(contact_num > 0, 1, 0)
        contacts_by_helix_frame_resi.append(contact_bin)
    
    #should be [n_helices x n_upperleaflet_aa x n_frames]
    #contacts_by_helix_frame_resi = np.stack(contacts_by_helix_frame_resi)
    #print(contacts_by_helix_frame_resi.shape) 

    #find frames/FA where the FA touches both neighboring helices

    contacts_both_neighbors = []
    for i in range(6):
        cbn = np.multiply(contacts_by_helix_frame_resi[i], contacts_by_helix_frame_resi[(i+1)%6])
        #print(contacts_twoadjacent.shape)
        #cbn = np.where(contacts_by_helix_frame_resi[i] == 1 and contacts_by_helix_frame_resi[(i+1)%6] == 1, 1, 0)
        #print(cbn.shape)
        contacts_both_neighbors.append(cbn)
    contacts_both_neighbors = np.stack(contacts_both_neighbors)

    contacts_and_planecrossing = np.multiply(contacts_both_neighbors, planes[0:6])

    return contacts_both_neighbors, contacts_and_planecrossing


identify_bound_states(_identify_bound_states)
