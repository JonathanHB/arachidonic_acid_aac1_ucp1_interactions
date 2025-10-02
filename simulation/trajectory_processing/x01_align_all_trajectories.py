import os

#TODO put all the output of this in its own folder for better organization and easier disposal 

wynton_path = "/media/X01Raid01/Data_Backup/shared/UCP1_simulations/long_nonadaptive_simulations/long-aac1-ucp1-simulations-wynton/"
degrabo_path = "/media/X01Raid01/Data_Backup/shared/UCP1_simulations/long_nonadaptive_simulations/long-aac1-ucp1-simulations-degrabo/"

output_path = "/media/X01Raid01/Data_Backup/home/jborowsky/long-sims/long-aac1-ucp1-processing/centered_trajectories"

proteins = ["aac1", "ucp1"]

for p in proteins:

    #wynton data processing
    folders = os.listdir(wynton_path+p)
    runs = [f for f in folders if f[0:3] == "run"]

    for r in runs:
        path = wynton_path+p+"/"+r
        os.chdir(path)

        #figure out how many complete sets of 100 trajectory segments there are
        xtcfiles = [f for f in os.listdir() if (f[0:14] == "traj_comp.part" and f[-4:] == ".xtc")]
        xtcinds = [int(xf[14:18]) for xf in xtcfiles]
        
        if max(xtcinds)%100 == 99:
            c_max = max(xtcinds)//100
        else:
            c_max = max(xtcinds)//100 #-1 #the minus 1 makes the script not process trajectories until a full 100 (or 99) are present

        #concatenate and align each set of 100 trajectory segments
        for c in range(c_max+1):
            
            #concatenate
            if c == 0:
                catpath = f"{output_path}/wynton/{p}/{r}/trj-{c}01-{c}99"
            elif c == c_max:
                catpath = f"{output_path}/wynton/{p}/{r}/trj-{c}00-{max(xtcinds)}"
            else:
                catpath = f"{output_path}/wynton/{p}/{r}/trj-{c}00-{c}99"

            if not os.path.exists(catpath+".xtc") and not os.path.exists(catpath+"-centered.xtc"):
                g1 = os.system(f"gmx trjcat -f traj_comp.part0{c}*.xtc -o {catpath}.xtc")
                
            #center and make molecules whole
            if not os.path.exists(catpath+"-centered.xtc"):
                g1 = os.system(f"echo 1 0 | gmx trjconv -f {catpath}.xtc -s ../input/{p}_500ns.tpr -o {catpath}-centered.xtc -pbc mol -center")
                #os.system(f"rm {catpath}.xtc")

            #align and recenter
            if not os.path.exists(catpath+"-aligned-xy-s10.xtc"):
                continue
                g1 = os.system(f"echo 1 1 0 | gmx trjconv -f {catpath}-centered.xtc -s ../input/{p}_500ns.tpr -o {catpath}-aligned-xy-s10.xtc -fit rotxy+transxy -center -skip 10")
                #os.system(f"rm {catpath}-centered.xtc")

            #switch to pbc atom
            if not os.path.exists(catpath+"-aligned-xy-s10-atom.xtc"):
                continue
                g1 = os.system(f"echo 1 0 | gmx trjconv -f {catpath}-aligned-xy-s10.xtc -s ../input/{p}_500ns.tpr -o {catpath}-aligned-xy-s10-atom.xtc -center -pbc atom")
                


    #degrabo data processing
    folders = os.listdir(degrabo_path+p)
    runs = [f for f in folders if f[0:3] == "run"]
    
    for r in runs:
        path = degrabo_path+p+"/"+r
        os.chdir(path)
        
        #figure out how many complete sets of 100 trajectory segments there are
        xtcfiles = [f for f in os.listdir() if (f[0:14] == "traj_comp.part" and f[-4:] == ".xtc")]

        for xtcfile in xtcfiles:

            fnbase = f"{output_path}/degrabo/{p}/{r}/{xtcfile[0:18]}"
            #center and make molecules whole
            if not os.path.exists(fnbase+"-centered.xtc"):
                g1 = os.system(f"echo 1 0 | gmx trjconv -f {xtcfile} -s ../input/input_200ns.tpr -o {fnbase}-centered.xtc -pbc mol -center")

            #align and recenter
            if not os.path.exists(fnbase+"-aligned.xtc"):
                continue
                g1 = os.system(f"echo 1 1 0 | gmx trjconv -f {fnbase}-centered.xtc -s ../input/input_200ns.tpr -o {fnbase}-aligned-xy.xtc -fit rotxy+transxy -center")
                os.system(f"rm {fnbase}-centered.xtc")

# gmx trjcat -f traj_comp.part0*.xtc -o trj-001-110.xtc
# gmx trjconv -f trj-001-110.xtc -s ../input/aac1_500ns.tpr -o trj-001-110-centered.xtc -pbc mol -center
# gmx trjconv -f trj-001-110-centered.xtc -s ../input/aac1_500ns.tpr -o trj-001-110-aligned.xtc -fit rot+trans -center

#note that -pbc mol and -fit rot+trans don't go together in one command
