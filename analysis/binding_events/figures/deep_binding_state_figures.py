import pymol
from pymol import cmd
#import cmd.util
#the linter is wrong about the imports here

import os

resns = {"aa":"ARA", "popc":"POP"}
views = {"ucp1":( \
            0.420852005,   -0.896865129,   -0.136042103,\
            -0.899247408,   -0.432184339,    0.067328125,\
            -0.119178399,    0.094002955,   -0.988402486,\
            0.001598268,    0.000872221, -139.467361450,\
            50.793403625,   48.103023529,   43.123470306,\
            66.959152222,  212.183853149,  -20.000000000 ), \
         "aac1": ( \
            0.420852005,   -0.896865129,   -0.136042103,\
            -0.899247408,   -0.432184339,    0.067328125,\
            -0.119178399,    0.094002955,   -0.988402486,\
            0.001598268,    0.000872221, -139.467361450,\
            50.793403625,   48.103023529,   43.123470306,\
            66.959152222,  212.183853149,  -20.000000000 )}


molecule = "aa"
proteins = ["ucp1", "aac1"]

ref_path = "/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/wynton/ucp1/input/seg_0035.gro"
gro_path = "/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/long-aac1-ucp1-processing/wynton"
upperpath = "/home/jonathan/Documents/grabelab/aac1-ucp1/long-aac1-ucp1/closest_approaches/figure_frames"
outpath = "/home/jonathan/Documents/grabelab/aac1-ucp1/paper-01/independent-binding-events/figures/pngs"

if molecule == "aa":
    files = {"aac1":[
                ["cavity", "aac1-wynton-run04-initseg200-ARAN-rseq314-upper-frame5806-0.1nm"], \
                ["h1h2", "aac1-wynton-run04-initseg200-ARAN-rseq312-upper-frame7801-0.91nm"], \
                ["h3h4", "aac1-wynton-run02-initseg001-ARAN-rseq314-upper-frame3850-1.11nm"], \
                ["h5h6", "aac1-wynton-run02-initseg200-ARAN-rseq316-upper-frame2303-1.03nm"], \
                ["h5h6", "aac1-wynton-run03-initseg400-ARAN-rseq311-upper-frame7590-0.12nm"]],\
            "ucp1":[ \
                ["h5h6", "ucp1-degrabo-run01-initseg005-ARAN-rseq533-upper-frame594-0.49nm"], \
                ["h3h4", "ucp1-degrabo-run01-initseg014-ARAN-rseq527-upper-frame753-0.14nm"]]}


elif molecule == "popc":
    files = {"aac1":[
                ["h1h2", "aac1-wynton-run01-initseg001-POPC-rseq425-upper-frame6189-0.8nm"], \
                ["h6h1", "aac1-degrabo-run02-initseg009-POPC-rseq448-upper-frame1534-0.4nm"], \
                ["h5h6", "aac1-wynton-run04-initseg200-POPC-rseq454-upper-frame8718-0.4nm"]],\
            "ucp1":[ \
                ["h5h6", "ucp1-degrabo-run01-initseg003-POPC-rseq411-upper-frame18-0.81nm"], \
                ["h3h4", "ucp1-degrabo-run01-initseg014-POPC-rseq423-upper-frame87-1.11nm"]]}
                

for protein in proteins:
    for f in files[protein]:

        site = f[0]
        fn = f[1]

        cmd.delete("all")

        #load and align
        cmd.load(ref_path, "ref")
        cmd.load(f"{gro_path}/{protein}/input/seg_0035.gro")
        cmd.cealign("ref and poly and name CA", "seg_0035") #make alignment consistent between aac1 and ucp1


        cmd.load(f"{upperpath}/{molecule}/{protein}/{site}/{fn}.pdb", "frame")

        cmd.cealign("seg_0035 and poly and name CA", "frame")

        #display
        cmd.hide("everything")

        cmd.dss()
        cmd.show("cartoon", "frame")

        resi = int(fn.split("/")[-1].split("-")[5][4:])
        #print(resi)

        cmd.show("spheres", f"frame and resn {resns[molecule]} and resi {resi} and not elem H")
        
        threshold = 5

        cmd.show("sticks", f"byres((frame and poly) within {threshold} of (frame and resn {resns[molecule]} and resi {resi}))")

        cmd.hide("sticks", "elem H")

        #coloring
        util.cbag()
        util.cbac(f"byres((frame and poly) within {threshold} of (frame and resn {resns[molecule]} and resi {resi}))")
        util.cbay(f"frame and resn {resns[molecule]} and resi {resi}")

        cmd.set_view((\
     0.420852005,   -0.896865129,   -0.136042103,\
    -0.899247408,   -0.432184339,    0.067328125,\
    -0.119178399,    0.094002955,   -0.988402486,\
     0.001598268,    0.000872221, -181.657089233,\
    50.793403625,   48.103023529,   43.123470306,\
   109.148910522,  254.373580933,  -20.000000000 ))
        

        cmd.png(f"{outpath}/{fn}-{site}-yellow.png", width=2400, height=1800, ray=True)
        #sys.exit(0)

