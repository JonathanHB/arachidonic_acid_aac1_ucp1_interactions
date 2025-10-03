#Jonathan Borowsky
#Grabe lab
#8/6/25

#run this script with: source /home/jonathan/Documents/grabelab/aac1-ucp1/paper-01/processing-and-visualization/load_degrabo_trajectories_v1.tcl

mol new /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/degrabo/aac1/input/seg_0035.gro type gro

set trj_upper /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/degrabo/aac1/run01/
#note these lists are not comma-delimited
set trajectories {traj_comp.part0001-aligned.xtc traj_comp.part0002-aligned.xtc traj_comp.part0003-aligned.xtc traj_comp.part0004-aligned.xtc traj_comp.part0005-aligned.xtc traj_comp.part0006-aligned.xtc traj_comp.part0007-aligned.xtc traj_comp.part0008-aligned.xtc traj_comp.part0009-aligned.xtc traj_comp.part0010-aligned.xtc traj_comp.part0011-aligned.xtc}

for {set i 0} {$i < 11} {incr i} {
 
    mol addfile ${trj_upper}[lindex $trajectories $i] type xtc step 20
    echo ${trj_upper}[lindex $trajectories $i]
}

set j 0

mol addrep $j
mol modselect 1 $j all and not name \"C.*\" and not name \"H.*\"
mol modstyle 1 $j Lines
mol modcolor 1 $j name

mol addrep $j
mol modselect 2 $j protein
mol modstyle 2 $j Ribbons
mol modcolor 2 $j ColorID 7

mol addrep $j
mol modselect 3 $j resname ARAN and not name \"H.*\"
mol modstyle 3 $j Licorice 0.300000 12.000000 12.000000
mol modcolor 3 $j resID