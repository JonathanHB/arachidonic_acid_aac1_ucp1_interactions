import os

toppath_upper = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories"
inputpath     = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories_links"

servers = ["wynton", "degrabo"]
proteins = ["aac1", "ucp1"]

def compute_observable_x01(calc_obs):

    #loop over servers, proteins, and parallel runs
    for server in servers:
        for protein in proteins:
            toppath = f"{toppath_upper}/{server}/{protein}/input/seg_0035.gro"
            for run in range(1,5):
                #not all servers have the same number of parallel simulations
                if run > 2 and server == "degrabo": continue

                for seg in range(1,20):
                    segstr = str(seg).zfill(2)

                    trjpath = f"{inputpath}/{server}/{protein}/run0{run}/seg{segstr}.xtc"
                    if not os.path.exists(trjpath): break
                    
                    #upper level directory information is supplied by the calc_obs function
                    outpath = f"{server}/{protein}/run0{run}/{protein}-{server}-run0{run}-seg{segstr}"
                    calc_obs(trjpath, toppath, protein, outpath)
                    print(f"saved output for {[server, protein, run, seg]}")
                
