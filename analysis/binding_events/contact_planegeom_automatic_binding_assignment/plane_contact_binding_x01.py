import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import mdtraj as md


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

    resqueries = []
    for rt in terminal_residues[protein]:
        helix_resseqs = [i for i in range(rt[0], rt[1]+1)]
        #print("color red, resi " + "+".join([str(i) for i in helix_resseqs]))
        resqueries.append(" or ".join([f"resSeq {i}" for i in helix_resseqs]))
    
    return resqueries


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
        rq = helix_resqueries(protein)

        ref_trj = md.load(f"{toppath_upper}/{servers[0]}/{protein}/input/seg_0035.gro")
        ref_inds = ref_trj.top.select("protein and not element H")
        trj = md.load(f"{toppath_upper}/{servers[0]}/{protein}/input/seg_0035.gro", atom_indices = ref_inds)

        cumulative_atoms = 0
        helix_inds = []
        for rqi in rq:
            sele = trj.top.select(rqi)
            helix_inds.append([i+cumulative_atoms for i in range(len(sele))])
            cumulative_atoms += len(sele)

        #-----------------------------------------------------------------------------------------------------------------------------

        rn = 1
        for server in servers:
            for run in range(1,5):
                
                print(protein, run, server)

                planes_all = []
                contacts_all = []
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
                        #contacts.transpose(2,0,1)
                        #planes: [n_helices+1 x n_upperleaflet_aa x n_frames] 
                        planes = np.load(file_planes).astype("uint8")
                        
                        cbn, cap = ibs(contacts, planes, helix_inds)
                        #sys.exit(0)
                        planes_all.append(planes)
                        contacts_all.append(cbn)
                        contacts_planes_all.append(cap)
                        
                # print(len(planes_all))
                planes_all = np.concatenate(planes_all, axis=2)
                contacts_all = np.concatenate(contacts_all, axis=2)
                contacts_planes_all = np.concatenate(contacts_planes_all, axis=2)
                
                print(planes_all.shape)
                print(contacts_all.shape)
                print(contacts_planes_all.shape)

                colors = ["purple", "blue", "cyan", "green", "orange", "red", "magenta", "grey"]
                legend = ["1", "2", "3", "4", "5", "6", "7", "8"]

                plt.figure(dpi=600)
                for i in range(6):
                    #plt.title(f"h{i+1}-h{(i+1)%6+1}")
                    for j in range(8):
                        #plt.plot(planes_all[i][j], color = "blue")
                        #plt.plot(contacts_all[i][j], color = "red")
                        plt.scatter([t*step for t in range(len(contacts_planes_all[i][j]))], (i+1)*contacts_planes_all[i][j]+j/32-1/8, color = colors[j], s=0.5, marker="o")

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
