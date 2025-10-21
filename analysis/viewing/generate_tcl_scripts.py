import os
import sys

current_path = "/home/jonathan/Documents/grabelab/aac1-ucp1/arachidonic_acid_aac1_ucp1_interactions/analysis/viewing/"
trj_upper = "/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/trajectories_for_viewing/aligned_trajectories_s20/"

tcl_contents = """
set j 0

#mol addrep $j
mol modselect 0 $j all and not name \\"C.*\\" and not name \\"H.*\\"
mol modstyle 0 $j Lines
mol modcolor 0 $j name

mol addrep $j
mol modselect 1 $j protein
mol modstyle 1 $j Ribbons
mol modcolor 1 $j ColorID 7

mol addrep $j
mol modselect 2 $j resname ARAN and not name \\"H.*\\"
mol modstyle 2 $j Licorice 0.300000 12.000000 12.000000
mol modcolor 2 $j Name #resID"""


for protein in ["aac1", "ucp1"]:
    ri = 1
    for server in ["wynton", "degrabo"]:
        for run in range(1,5):
            if server == "degrabo" and run > 2:
                continue

            runstr = str(run).zfill(2)
            trjpath = f"{trj_upper}{server}/{protein}/run{runstr}/"
            tcl_fn = f"load-{protein}-run{ri}.tcl"

            trjfiles = []

            nsegs = 0
            for seg in range(20):

                segstr = str(seg).zfill(2)
                segfn = f"{protein}-{server}-run{runstr}-seg{segstr}-aligned-s20.xtc"
                # print(trjpath+segfn)
                if os.path.exists(trjpath+segfn):
                    trjfiles.append(segfn)
                    nsegs += 1

            trjlist = " ".join(trjfiles)


            file_contents = [
                f"# run this file in VMD TK Console with: source {current_path+tcl_fn}\n",
                f"mol new {trj_upper}{protein}_input/seg_0035.gro type gro",
                f"set trj_upper {trjpath}",
                f"set trajectories {{{trjlist}}}",
                f"\nfor {{set i 0}} {{$i < {nsegs}}} {{incr i}} {{",
                "   mol addfile ${trj_upper}[lindex $trajectories $i] type xtc step 1",
                "   echo ${trj_upper}[lindex $trajectories $i]",
                "}",
                tcl_contents
            ]

            tcl_fn = f"load-{protein}-run{ri}.tcl"
            with open(tcl_fn, "w") as f:
                f.writelines(ln+"\n" for ln in file_contents)

            f.close()

            ri += 1