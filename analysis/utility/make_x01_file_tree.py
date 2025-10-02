import os

for server in ["wynton", "degrabo"]:
    os.mkdir(server)
    for protein in ["aac1", "ucp1"]:
        os.mkdir(f"{server}/{protein}")
        for run in range(1,5):
            os.mkdir(f"{server}/{protein}/run0{run}")
            
