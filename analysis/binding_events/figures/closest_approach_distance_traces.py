import numpy as np
# import mdtraj as md
#from protein_ligand_contacts import get_contacts
#from protein_ligand_contacts import save_pdb_bfactors
# import time
# from enspara.info_theory.exposons import exposons_from_sasas
import os
import sys

import matplotlib.pyplot as plt

# sys.path.insert(1, f'{os.getcwd()}/../utility')

# import load_all_trjs
# import load_anatale_trjs
# import itertools
# import count_lipids

# import compute_observables_on_trjs_x01_v2

tracked_aa_upper = [["ucp1", "wynton", 1, 530, "h56"],
                    ["ucp1", "wynton", 3, 533, "h56"],
                    ["ucp1", "wynton", 4, 528, "h56"],
                    ["ucp1", "wynton", 4, 529, "h56"],
                    ["ucp1", "degrabo", 1, 527, "h34"],
                    ["ucp1", "degrabo", 1, 533, "h56"],

                    ["aac1", "degrabo", 1, 309, "h12"],
                    ["aac1", "degrabo", 1, 311, "h12"],
                    ["aac1", "wynton", 3, 313, "h12"],
                    ["aac1", "wynton", 3, 315, "h12"],
                    ["aac1", "wynton", 4, 312, "h12"],

                    ["aac1", "wynton", 2, 314, "h34"],
                    ["aac1", "wynton", 3, 314, "h34"],

                    ["aac1", "wynton", 1, 311, "h56"],
                    ["aac1", "wynton", 2, 315, "h56"],
                    ["aac1", "wynton", 2, 316, "h56"],
                    ["aac1", "wynton", 3, 311, "h56"],

                    ["aac1", "degrabo", 1, 316, "cavity"],
                    ["aac1", "degrabo", 2, 310, "cavity"],
                    ["aac1", "wynton", 4, 314, "cavity"]
                    ]


def compile_distance_data(residue, leaflet):
    
    stride = 1
    window = 100

    inputpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/closest_approaches"

    servers = ["wynton", "degrabo"]
    proteins = ["aac1", "ucp1"]

    #loop over servers, proteins, and parallel runs
    for protein in proteins:

        dist_traces = {"h12":[], "h34":[], "h56":[], "cavity":[]}
        dist_legends = {"h12":[], "h34":[], "h56":[], "cavity":[]}

        run_num = 0

        pml_out_lines = ["delete all\n", f"load /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/wynton/{protein}/input/seg_0035.gro\n"]
        for server in servers:
            #toppath = f"{inputpath}/{server}/{protein}/input/seg_0035.gro"
            for run in range(1,5):
                #print(run)
                #segments are named differently on different servers
                if server == "wynton":
                    segs = ["001"]+[str(i) for i in range(100, 1000, 100)]
                elif server == "degrabo":
                    segs = [str(i).zfill(3) for i in range(1,20)]
                else:
                    segs = [] #ensure variable is defined

                #step is in microseconds
                if server == "wynton" and (run == 3 or run == 4):
                    step = 1/10000
                else:
                    step = 1/5000


                fnbase = f"{inputpath}/{server}/{protein}/run0{run}/{protein}-{server}-run0{run}"
                if not os.path.exists(f"{fnbase}-initseg{segs[0]}-{residue}-{leaflet}-resseqs.npy"):
                    print(f"no file found for {server} {protein} run{run}")
                    continue
                
                #print(server, protein, run, leaflet, residue)
                resseqs = list(np.load(f"{fnbase}-initseg{segs[0]}-{residue}-{leaflet}-resseqs.npy"))
                #print(resseqs)
                #figure out which trajectory segment has the closest approach
                #min_dists_all = []
                #argmin_dists_all = []
                dists_all = []

                for seg in segs:
                    spath = f"{fnbase}-initseg{seg}-{residue}-{leaflet}.npy"
                    if os.path.exists(spath):
                        dists = np.load(spath)
                        dists_all.append(dists)
                        #min_dists_all.append(np.min(dists, axis = 1))
                        #argmin_dists_all.append(np.argmin(dists, axis = 1))
                        
               	#for d in dists_all:
                #    print(d.shape) 
                
                dists_all = np.concatenate(dists_all, axis = 1)
                
                for taa in tracked_aa_upper:
                    if taa[0] == protein and taa[1] == server and taa[2] == run:
                        dist_trj = dists_all[resseqs.index(str(taa[3]))]
                        dist_smoothed = [np.mean(dist_trj[i:i+window]) for i in range(len(dist_trj)-window)]
                        dist_traces[taa[4]].append([[i*step for i in range(dists_all.shape[1]-window)], dist_smoothed, dist_trj[:-window]]) #, linewidth = 0.5)
                        #dist_traces.append(dists_all[resseqs.index(taa[3])])
                        dist_legends[taa[4]].append(f"run {1+run_num}, AA {1+resseqs.index(str(taa[3]))}")

                run_num += 1

        for site in ["h12","h34","h56","cavity"]:
            
            if len(dist_traces[site]) == 0:
                continue
	
            colors = ["blue", "orange", "green", "red", "purple", "yellow"]

            plt.figure()

            for ci, data in enumerate(dist_traces[site]):
                plt.plot(data[0], data[2], linewidth = 0.5, alpha = 0.25, color = colors[ci], label="_nolegend_")
                plt.plot(data[0], data[1], color = colors[ci], label = f"line{ci}")

            #plt.title(f"{protein} {site}")
            plt.xlabel("time (microseconds)")
            plt.ylabel(f"AA COO- distance to protein center of mass (nm)")
            plt.legend(dist_legends[site])
            plt.ylim(0,3)
            plt.xlim(0,7.5)
            svg_path = f"/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/closest_approaches/closest_approaches_{protein}_{residue}_{leaflet}_{site}.svg"
            plt.savefig(svg_path, format='svg')
            #png_path = f"/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/closest_approaches/closest_approaches_{protein}_{residue}_{leaflet}_{site}.png"
            #plt.savefig(png_path, dpi = 600)
            plt.clf()

                #print(dists_all.shape)
                #sys.exit(0)



                # argmin_dists_all = np.stack(argmin_dists_all)

                # argmin_segs = np.argmin(min_dists_all, axis = 0)
                # print(argmin_segs)
                #assemble pml script to load the closest approach of each residue in the current run
                #pml_out_lines = ["delete all\n"] #, f"load {top_path}\n"]

                # for ai, ams in enumerate(argmin_segs):
                    
                #     dist_threshold = 1.3
                #     if min_dists_all[ams][ai] > dist_threshold:
                #         continue
                #     #print(ai)
                #     #print(ams)
                #     #print(segs[ams])
                #     #print(resseqs[ai])
                #     #print(argmin_dists_all[ams][ai])
                #     #print(round(min_dists_all[ams][ai], 2))
                #     #print("----------------------------------------------")
                #     dist_str = str(round(min_dists_all[ams][ai], 2))
                    
                #     outputpath = "/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/closest_approaches"
                #     fnout = f"{outputpath}/{server}/{protein}/run0{run}/{protein}-{server}-run0{run}"
                #     savefn = f"{fnout}-initseg{segs[ams]}-{residue}-rseq{resseqs[ai]}-{leaflet}-frame{argmin_dists_all[ams][ai]*stride}-{dist_str}nm.pdb"
                #     protfn = savefn.split("/")[-1][:-4]

                #     pml_out_lines.append(f"load {savefn}\n")
                #     pml_out_lines.append(f"cealign seg_0035 and poly and name CA, {protfn}\n")
                    
                #     pml_out_lines.append(f"hide everything, {protfn}\n")
                #     pml_out_lines.append(f"show sticks, {protfn} and resi {resseqs[ai]}\n")

#         #add graphics commands and save script
#         pml_out_lines.append(f"show cart\n")
#         pml_out_lines.append(f"show lines, resn ARA\n")
#         pml_out_lines.append(f"util.cbag\n")
#         pml_out_lines.append(f"util.cbac resn POPC\n")
#         pml_out_lines.append(f"util.cbam resn ARA\n")

#         pml_path = f"/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/closest_approaches/load_closest_approaches_{protein}_{residue}_{leaflet}.pml"

#         if os.path.exists(pml_path): #not redundant with "w" if the original file is longer than the new one
#             os.remove(pml_path)
#         with open(pml_path, "w") as f:
#             f.writelines(pml_out_lines)
#         f.close()


# for leaflet in ["upper", "lower"]:
#     for residue in ["ARAN", "POPC"]:
compile_distance_data("POPC", "upper")
