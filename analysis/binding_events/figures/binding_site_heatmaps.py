import pymol
from pymol import cmd

import os


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

        resqueries.append("+".join([str(i) for i in helix_resseqs]))

    return resqueries

upperpath = f"/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/contact_freq_heatmaps"

proteins = ["ucp1", "aac1"]

for protein in proteins:

    rq = helix_resqueries(protein)

    for si in range(6):
        fpath = f"{upperpath}/bound-contacts-{protein}-site-{si}.pdb"
        if os.path.exists(fpath):

            cmd.delete("all")

            cmd.load(fpath, f"{protein}-site-{si}")
            cmd.hide("everything", f"{protein}-site-{si}")
            cmd.show("cart", f"{protein}-site-{si}")
            cmd.spectrum("b", "blue_red")
            cmd.show("sticks", f"resi {rq[si]} and poly and not elem H")
            cmd.show("sticks", f"resi {rq[((si+1)%6)]} and poly and not elem H")

            cmd.color("grey", f"not resi {''.join(rq)}")

            if protein == "aac1":
                if si == 5:
                    cmd.set_view((\
                        0.198175028,   -0.051406000,   -0.978818119,\
                        0.979623020,   -0.022814510,    0.199535340,\
                        -0.032586936,   -0.998409927,    0.045835838,\
                        -0.000110248,    0.000507869, -164.618713379,\
                        42.775024414,   54.231544495,   46.452850342,\
                        139.773864746,  189.477798462,  -20.000000000 ))
                elif si == 4:
                    cmd.set_view((\
                        -0.932600617,   -0.064751804,   -0.355044246,\
                        0.353186280,    0.038548052,   -0.934753180,\
                        0.074212536,   -0.997143626,   -0.013081808,\
                        -0.000014657,    0.001085736, -173.493804932,\
                        46.247146606,   46.133419037,   46.630638123,\
                        152.120605469,  194.897598267,  -20.000000000 ))
                elif si == 2:
                    cmd.set_view((\
                        -0.136997581,   -0.023314206,    0.990293086,\
                        -0.982807159,   -0.121699497,   -0.138826832,\
                        0.123753332,   -0.992279232,   -0.006239942,\
                        0.000259329,    0.001263648, -182.883728027,\
                        53.974460602,   49.783386230,   45.033416748,\
                        160.280029297,  205.480621338,  -20.000000000 ))
                    # cmd.set_view((\
                    #     -0.136997581,   -0.023314206,    0.990293086,\
                    #     -0.982807159,   -0.121699497,   -0.138826832,\
                    #     0.123753332,   -0.992279232,   -0.006239942,\
                    #     0.000259329,    0.001263648, -182.883728027,\
                    #     53.974460602,   49.783386230,   45.033416748,\
                    #     111.927040100,  253.833633423,  -20.000000000 ))
                elif si == 0:
                    cmd.set_view((\
                        0.939353943,   -0.063957915,   -0.336872041,\
                        0.329026967,   -0.024412701,    0.942972302,\
                        -0.049055371,   -0.998732090,   -0.008410829,\
                        0.653967619,   -1.351881266, -213.204727173,\
                        64.210937500,    7.937337875,   48.915527344,\
                        123.683120728,  186.953109741,  -20.000000000 ))
            
            elif protein == "ucp1":
                if si == 4:
                    cmd.set_view((\
                        -0.952760696,   -0.245877430,    0.178281546,\
                        -0.222454607,    0.165315479,   -0.960818648,\
                        0.206768230,   -0.955085754,   -0.212203115,\
                        0.000039084,    0.001076169, -168.115539551,\
                        54.714576721,   37.325725555,   44.832897186,\
                        148.818817139,  187.448272705,  -20.000000000 ))
                    # cmd.set_view((\
                    #     -0.952760696,   -0.245877430,    0.178281546,\
                    #     -0.222454607,    0.165315479,   -0.960818648,\
                    #     0.206768230,   -0.955085754,   -0.212203115,\
                    #     -0.000113636,    0.001297221, -174.137054443,\
                    #     53.588649750,   48.216567993,   47.844478607,\
                    #     153.004928589,  195.305694580,  -20.000000000 ))
                    # cmd.set_view((\
                    #     -0.965648770,   -0.245089814,    0.086301126,\
                    #     -0.110812329,    0.088013269,   -0.989932835,\
                    #     0.235024616,   -0.965482354,   -0.112150207,\
                    #     -0.000113636,    0.001297221, -174.137054443,\
                    #     53.588649750,   48.216567993,   47.844478607,\
                    #     153.004928589,  195.305694580,  -20.000000000 ))
                elif si == 2:
                    cmd.show("sticks", f"resi 95 and poly and not elem H")

                    cmd.set_view((\
                        0.661948681,   -0.305526853,    0.684450805,\
                        -0.716476381,    0.010317598,    0.697529972,\
                        -0.220175281,   -0.952114463,   -0.212068528,\
                        0.000265977,    0.000888564, -153.720764160,\
                        60.762977600,   57.539386749,   40.991134644,\
                        122.920829773,  184.532653809,  -20.000000000 ))
                    
            #this print statement may be important for getting the file to save due to some bug where pymol doesn't wait correctly
            print(f"{upperpath}/bound-contacts-{protein}-site-{si}-heatmap.png")
            cmd.png(f"{upperpath}/bound-contacts-{protein}-site-{si}-heatmap.png", ray = 1)


#color bar
cmd.delete("all")
cmd.ramp_new("colorbar", "none", [0,1], ["blue", "red"])
cmd.png(f"{upperpath}/bound-contacts-colorbar.png", ray = 1)