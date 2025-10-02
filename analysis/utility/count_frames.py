import os
import mdtraj as md

inputpath = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories"

def observable(calc_obs, server, protein, run, seg, trjpath, toppath):
    #outpath = f"{server}/{protein}/run0{run}/{protein}-{server}-run0{run}-n-frames.npy"
    trj = md.load(trjpath, top=toppath, atom_indices=[0])
    return trj.n_frames#, outpath
    #calc_obs(trjpath, toppath, protein, outpath)
    #print(f"saved output for {[server, protein, run, seg]}")

servers = ["wynton", "degrabo"]
proteins = ["aac1", "ucp1"]

def compute_observable_x01(calc_obs):

    data_lines = []

    #loop over servers, proteins, and parallel runs
    for server in servers:
        print(server)
        for protein in proteins:
            print(f"  {protein}")
            toppath = f"{inputpath}/{server}/{protein}/input/seg_0035.gro"
            for run in range(1,5):
                print(f"    {run}")
                #check if run actually exists since not all servers have the same number of parallel simulations
                fpath = f"{inputpath}/{server}/{protein}/run0{run}"
                if os.path.exists(fpath):
                    
                    total_frames = 0
                    #segments are named differently on different servers
                    if server == "wynton":
                        continue
                        #find all the trjcatted segments in order

                        xtcfiles = [f for f in os.listdir(fpath) if f[-13:] == "-centered.xtc"]
                        
                        for seg in ["001"]+[str(i) for i in range(100, 1000, 100)]:
                            xtcfile = [x for x in xtcfiles if x[4:7] == seg]
                            
                            if len(xtcfile) == 1:
                                total_frames += observable(calc_obs, server, protein, run, seg, f"{fpath}/{xtcfile[0]}", toppath)
                                # outpath = f"{server}/{protein}/run0{run}/{protein}-{server}-run0{run}-initseg{str(seg).zfill(3)}"
                                # calc_obs(xtcfile[0], toppath, protein, outpath)
                                # print(f"saved output for {[server, protein, run, seg]}")
                                #total_frames += frames

                    
                    if server == "degrabo":
                        
                        #the segments here are more standardized than on wynton so it's easier to loop over them
                        for seg in range(20):
                            trjpath = f"{fpath}/traj_comp.part{str(seg).zfill(4)}-centered.xtc"
                            if os.path.exists(trjpath):
                                total_frames += observable(calc_obs, server, protein, run, seg, trjpath, toppath)
                                #total_frames += frames


                    data_lines.append(" ".join([str(s) for s in [protein, server, run, total_frames]]))

    print(data_lines)
    out_path = f"/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/frame_counts_degrabo.txt"

    if os.path.exists(out_path): #not redundant with "w" if the original file is longer than the new one
        os.remove(out_path)
    with open(out_path, "w") as f:
        f.writelines(data_lines)
    f.close()

compute_observable_x01(observable)
