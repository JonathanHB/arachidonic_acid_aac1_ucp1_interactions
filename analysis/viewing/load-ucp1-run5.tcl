# run this file in VMD TK Console with: source /home/jonathan/Documents/grabelab/aac1-ucp1/arachidonic_acid_aac1_ucp1_interactions/analysis/viewing/load-ucp1-run5.tcl

mol new /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/trajectories_for_viewing/aligned_trajectories_s20/ucp1_input/seg_0035.gro type gro
set trj_upper /home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/trajectories_for_viewing/aligned_trajectories_s20/degrabo/ucp1/run01/
set trajectories {ucp1-degrabo-run01-seg01-aligned-s20.xtc ucp1-degrabo-run01-seg02-aligned-s20.xtc ucp1-degrabo-run01-seg03-aligned-s20.xtc ucp1-degrabo-run01-seg04-aligned-s20.xtc ucp1-degrabo-run01-seg05-aligned-s20.xtc ucp1-degrabo-run01-seg06-aligned-s20.xtc ucp1-degrabo-run01-seg07-aligned-s20.xtc ucp1-degrabo-run01-seg08-aligned-s20.xtc ucp1-degrabo-run01-seg09-aligned-s20.xtc ucp1-degrabo-run01-seg10-aligned-s20.xtc ucp1-degrabo-run01-seg11-aligned-s20.xtc ucp1-degrabo-run01-seg12-aligned-s20.xtc ucp1-degrabo-run01-seg13-aligned-s20.xtc ucp1-degrabo-run01-seg14-aligned-s20.xtc ucp1-degrabo-run01-seg15-aligned-s20.xtc}

for {set i 0} {$i < 15} {incr i} {
   mol addfile ${trj_upper}[lindex $trajectories $i] type xtc step 1
   echo ${trj_upper}[lindex $trajectories $i]
}

set j 0

#mol addrep $j
mol modselect 0 $j all and not name \"C.*\" and not name \"H.*\"
mol modstyle 0 $j Lines
mol modcolor 0 $j name

mol addrep $j
mol modselect 1 $j protein
mol modstyle 1 $j Ribbons
mol modcolor 1 $j ColorID 7

mol addrep $j
mol modselect 2 $j resname ARAN and not name \"H.*\"
mol modstyle 2 $j Licorice 0.300000 12.000000 12.000000
mol modcolor 2 $j Name #resID
