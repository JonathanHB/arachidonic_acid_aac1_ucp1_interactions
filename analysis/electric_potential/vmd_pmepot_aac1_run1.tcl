#Jonathan Borowsky
#Grabe lab
#9/20/25

#run this script with: source /home/jonathan/Documents/grabelab/aac1-ucp1/paper-01/electric_potential/092025/vmd_pmepot.tcl

mol new /home/jonathan/Documents/grabelab/aac1-ucp1/paper-01/electric_potential/092025/aac1/charge_addition/seg_0035.pqr type pqr
#mol addfile /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/wynton/aac1/input/seg_0035.gro type gro
#Note: pqr files yield strange bonds when loaded into vmd without a gro file to go with them.
#  This is normal (the pqr file contains no topology information of its own): https://dasher.wustl.edu/chem478/labs/lab-08/tutorial-vmd.pdf

set trj_upper /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/water_wires_2_noxyz/centered_trajectories_wire_s10/wynton/aac1/run01/
#note these lists are not comma-delimited
set trajectories {aac1-wynton-run01-seg01-wireframes-s10.xtc aac1-wynton-run01-seg02-wireframes-s10.xtc aac1-wynton-run01-seg03-wireframes-s10.xtc aac1-wynton-run01-seg04-wireframes-s10.xtc}

for {set i 0} {$i < 4} {incr i} {
    #without 'waitfor all' vmd will try to run pmepot on frames which have not yet loaded
    mol addfile ${trj_upper}[lindex $trajectories $i] type xtc step 1 waitfor all
    echo ${trj_upper}[lindex $trajectories $i]
}

#this seemingly must be run after loading trajectories and before pmepot to eliminate unwanted bonds; molid does not seem to work properly
#mol addfile /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/wynton/aac1/input/seg_0035.gro type gro molid 0

set sel_all [atomselect 0 "all"]

set outpath /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/electric_potential/test_aac1/

#start at 1 to avoid making an extra frame from the starting structure that doesn't correspond to a water wire
for {set i 1} {$i < 100} {incr i} {
    #set b [expr {$i}]
    pmepot -sel $sel_all -frames $i:$i -dxfile ${outpath}potential-wireframe-s10-$i.dx

    #outputs are too large
    #pmepot -sel $sel_all -frames $i:$i -dxfile ${outpath}potential-wireframe-s10-g120-120-110-$i.dx -grid 0.25

}

set j 0


#mol addrep $j
#mol modselect 1 $j all and not name \"C.*\" and not name \"H.*\"
#mol modstyle 1 $j Lines
#mol modcolor 1 $j name

#mol addrep $j
#mol modselect 2 $j protein
#mol modstyle 2 $j Ribbons
#mol modcolor 2 $j ColorID 7

#mol addrep $j
#mol modselect 3 $j resname ARAN and not name \"H.*\"
#mol modstyle 3 $j Licorice 0.300000 12.000000 12.000000
#mol modcolor 3 $j resID