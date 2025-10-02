#Jonathan Borowsky
#Grabe lab
#6/26/24s

# allow the script to execute pymol commands
from pymol import cmd
from pymol import util
# Import PyMOL's stored module.  This will allow us with a
# way to pull out the PyMOL data and modify it in our script.
from pymol import stored

cmd.delete("all")
cmd.reinitialize()

cmd.load("../seg_0035_ucp1_centered.gro")

cmd.create("protchain", "seg_0035_ucp1_centered and poly")
cmd.create("phosphates", "seg_0035_ucp1_centered and element P")

cmd.dss()

cmd.hide("everything", "object *")

cmd.set("specular", False)
cmd.color("lime")
cmd.bg_color("white")

cmd.show("surf", "object " + "protchain")
cmd.show("spheres", "object " + "phosphates")
cmd.color("orange", "phosphates")

cmd.set("surface_quality", 2)
cmd.set("ray_shadows", "off")
cmd.set("ray_interior_color", "gray50")

cmd.rotate("y",90)

cmd.set("orthoscopic", "on")

com = cmd.centerofmass("object protchain")

cmd.set_view((\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000001624,   -0.000001222, -281.098358154,\
    com[0],   com[1],   33.838760376,\
   272.977142334,  309.418609619,   20.000000000 ))


psi = 0
for i in range(-3,4):
    cmd.pseudoatom(f"pseudo_{psi}", pos=[com[0]+10*i,com[1],com[2]])
    cmd.color("black", f"pseudo_{psi}")
    cmd.show("nonbonded", f"pseudo_{psi}")

    psi+=1


cmd.png(f"/Users/jonathanborowsky/Documents/grabelab/aac1-ucp1/ucp1-cutaway_v1.png", width=2000, height=2000, ray=1)


#/Users/jonathanborowsky/Documents/grabelab/aac1-ucp1/paper-01/figures/scripts



# import sys
# sys.exit(0)

# object1 = "8hbvA"
# cmd.fetch(object1)
# cmd.set_name(object1, f"{object1}_ref")

# object1 = f"{object1}_ref"

# #reference residues
# bent_helix_m_ends = [42, 142, 241]
# bent_helix_c_ends = [12, 112, 211]
# straight_helix_m_ends = [77, 176, 270]
# straight_helix_c_ends = [94, 193, 287]
# #check that the same number of residues were included from each third of the protein
# print([i-j for i,j in zip(bent_helix_c_ends, bent_helix_m_ends)])
# print([i-j for i,j in zip(straight_helix_c_ends, straight_helix_m_ends)])

# helix_cterm_ends = bent_helix_m_ends+straight_helix_c_ends
# helix_nterm_ends = bent_helix_c_ends+straight_helix_m_ends

# # print(f"show cart, resi {' or resi '.join([str(i)+'-'+str(j) for i, j in zip(bent_helix_c_ends, bent_helix_m_ends)])}")
# # print(f"show cart, resi {' or resi '.join([str(i)+'-'+str(j) for i, j in zip(straight_helix_m_ends, straight_helix_c_ends)])}")

# #it turns out that the principal components of the well-ordered highly symmetric part of AAC1 are nearly degenerate,
# # so its symmetry axis does not lie along any of them

# #to get around this we calculate points on the symmetry axis
# # by averaging the coordinates of alpha carbons of pseudosymmetric trios of residues
# bent_helix_inds = [[k for k in range(i, j+1)] for i,j in zip(bent_helix_c_ends, bent_helix_m_ends)]
# straight_helix_inds = [[k for k in range(i, j + 1)] for i, j in zip(straight_helix_m_ends, straight_helix_c_ends)]

# all_reference_inds = []
# for i in bent_helix_inds+straight_helix_inds:
#     all_reference_inds += i

# #create pseudoatoms

# import numpy as np

# def make_symmetry_axis_pseudoatoms(atoms1, atoms2, atoms3, ps_index):
#     for i,j,k in zip(atoms1, atoms2, atoms3):

#         i_coords = cmd.get_coords(f"resi {i} and name CA", 1)[0]
#         j_coords = cmd.get_coords(f"resi {j} and name CA", 1)[0]
#         k_coords = cmd.get_coords(f"resi {k} and name CA", 1)[0]

#         axis_coords = np.mean(np.stack([i_coords, j_coords, k_coords]), axis=0)

#         cmd.pseudoatom(f"pseudo_{ps_index}", pos=list(axis_coords))

#         ps_index += 1

#     return ps_index

# #avoid overwriting pseudoatoms
# ps_ind = 0

# ps_ind = make_symmetry_axis_pseudoatoms(bent_helix_inds[0], bent_helix_inds[1], bent_helix_inds[2], ps_ind)
# cmd.color("red", "pseudo_*")
# ps_ind = make_symmetry_axis_pseudoatoms(straight_helix_inds[0], straight_helix_inds[1], straight_helix_inds[2], ps_ind)

# #cmd.hide("everything", "pseudo_*")


# #NOTE: including "object" in object selection prevents the bug where pymol doesn't recognize objects in selections
# #do not do this in fetch() commands

# object2 = "8hbvA"
# cmd.fetch(object2)
# cmd.align("object " + object2, f"object {object1} and resi {'+'.join([str(i) for i in all_reference_inds])}")

# #orient protein along its symmetry axis using pseudoatoms
# cmd.reset() #must be done after fetching since the fetched object shifts things
# cmd.center("pseudo_*")
# cmd.orient("pseudo_*")
# cmd.rotate("z", -90)
# cmd.rotate("y", 90)
# cmd.zoom("center", 35) #set for AAC1 M and C states in a 2000x2000 pixel image
# #note that this includes some sort of buffer

# cmd.dss()

# cmd.hide("everything", "object *")

# cmd.set("specular", False)
# cmd.color("lime")
# cmd.bg_color("white")

# # cmd.show("nonbonded", "object pseudo_*")
# # cmd.color("red", "object pseudo_*")

# util.cbay("resn GTP")
# cmd.show("sticks", "resn GTP")

# util.cbay("resn ATP")
# cmd.show("sticks", "resn ATP")

# util.cbay("resn DNF")
# cmd.show("sticks", "resn DNF")

# cutaway = 1
# # 0 = cartoon
# # 1 = cutaway
# # 2 = whole surface

# if cutaway == 0:
#     cmd.show("cartoon", "object " + object2)

#     #see https://pymolwiki.org/index.php/Cartoon_Helix_Settings
#     cmd.set("cartoon_fancy_helices", False)
#     cmd.set("cartoon_oval_length", 0.8)
#     cmd.set("cartoon_oval_width", 0.3)
#     cmd.set("cartoon_oval_quality", 100) #set to 1000 for final run?
#     cmd.set("cartoon_sampling", 20)


# if cutaway == 1:
#     cmd.show("surf", "object " + object2)

#     cmd.set("surface_quality", 2)
#     cmd.set("ray_shadows", "off")
#     cmd.set("ray_interior_color", "gray50")

#     manual_view = True

#     if not manual_view:
#         #clipping planes cannot be set differently for different objects
#         #these are ugly but keep the symmetry axis centered
#         if object2 == "8j1nA":
#             cmd.rotate("y", -90)
#             cmd.clip("near", -40)
#         elif object2 == "8hbwA":
#             cmd.rotate("y", -90)
#             cmd.clip("near", -40)
#             # cmd.rotate("y", 180)
#             # cmd.clip("near", -38.6) #0 deg and 39 or 270 deg and -38.6
#         elif object2 == "8hbvA":
#             cmd.rotate("y", -90)
#             cmd.clip("near", -40)

#     #These have been rotated to keep the projection of the symmetry axis vertical in the image.
#     # This guarantees nothing about whether it tilts into the plane of the screen.
#     else:
#         cmd.reset()

#         if object2 == "8j1nA":
#             cmd.set_view(( \
#                 -0.931455255, -0.247361913, 0.266828001, \
#                 -0.333892405, 0.872502267, -0.356713563, \
#                 -0.144571185, -0.421355098, -0.895296514, \
#                 0.000000000, 0.000000000, -201.273925781, \
#                 116.961257935, 121.269165039, 142.889236450, \
#                 199.273925781, 243.273941040, -20.000000000))

#             # cmd.set_view(( \
#             #     -0.263325721, 0.363038033, 0.893790007, \
#             #     -0.352128625, 0.826391935, -0.439403296, \
#             #     -0.898143709, -0.430435956, -0.089774221, \
#             #     0.000000000, 0.000000000, -201.273925781, \
#             #     116.961257935, 121.269165039, 142.889236450, \
#             #     199.273925781, 243.273941040, -20.000000000))

#         elif object2 == "8hbwA" or object2 == "8g8wA":
#             cmd.set_view(( \
#                 -0.139716908, -0.225761563, 0.964110076, \
#                 -0.346182168, 0.923352599, 0.166048065, \
#                 -0.927702844, -0.310556084, -0.207162619, \
#                 0.000001798, -0.000004835, -205.644180298, \
#                 116.410858154, 121.141448975, 143.195556641, \
#                 202.644180298, 249.813613892, -20.000000000))

#             # cmd.set_view(( \
#             #     -0.141680732, -0.224534348, 0.964110076, \
#             #     -0.338115394, 0.926337004, 0.166048065, \
#             #     -0.930376291, -0.302452624, -0.207162619, \
#             #     0.000001798, -0.000004835, -205.644180298, \
#             #     116.410858154, 121.141448975, 143.195556641, \
#             #     202.644180298, 249.813613892, -20.000000000))
#             # cmd.set_view(( \
#             #     -0.959719539, -0.241607472, 0.143397778, \
#             #     -0.269579858, 0.935656071, -0.227750868, \
#             #     -0.079147212, -0.257232755, -0.963102698, \
#             #     0.000019640, -0.000020198, -204.628921509, \
#             #     116.541992188, 120.900085449, 142.218185425, \
#             #     200.563308716, 248.894531250, -20.000000000))

#             # cmd.set_view(( \
#             #     -0.972342074, -0.184354469, 0.143397778, \
#             #     -0.213702679, 0.949977398, -0.227750868, \
#             #     -0.094240189, -0.252094686, -0.963102698, \
#             #     0.000019640, -0.000020198, -204.628921509, \
#             #     116.541992188, 120.900085449, 142.218185425, \
#             #     200.563308716, 248.894531250, -20.000000000))

#         elif object2 == "8hbvA":
#             cmd.set_view(( \
#                 -0.313199431, -0.430462301, -0.846525192, \
#                 0.231343210, 0.829937279, -0.507619739, \
#                 0.921074748, -0.354825407, -0.160351112, \
#                 0.000000000, 0.000000000, -219.326431274, \
#                 116.949813843, 120.683715820, 141.931060791, \
#                 217.326431274, 261.326416016, -20.000000000))

#             # cmd.set_view(( \
#             #     -0.294409275, -0.443525553, -0.846525192, \
#             #     0.195459738, 0.839113116, -0.507619739, \
#             #     0.935473561, -0.314910740, -0.160351112, \
#             #     0.000000000, 0.000000000, -219.326431274, \
#             #     116.949813843, 120.683715820, 141.931060791, \
#             #     217.326431274, 261.326416016, -20.000000000))

#     #cmd.ray()

# if cutaway == 2:
#     cmd.show("surf", "object " + object2)

#     cmd.set("surface_quality", 2)
#     cmd.set("ray_shadows", "off")
#     cmd.set("ray_interior_color", "gray50")

#     manual_view = False

#     if not manual_view:

#         pass

#         #clipping planes cannot be set differently for different objects
#         #these are ugly but keep the symmetry axis centered
#         # if object2 == "8j1nA":
#         #     cmd.rotate("y", -90)
#         #     cmd.clip("near", -40)
#         # elif object2 == "8hbwA":
#         #     cmd.rotate("y", -90)
#         #     cmd.clip("near", -40)
#         #     # cmd.rotate("y", 180)
#         #     # cmd.clip("near", -38.6) #0 deg and 39 or 270 deg and -38.6
#         # elif object2 == "8hbvA":
#         #     cmd.rotate("y", -90)
#         #     cmd.clip("near", -40)

#     #These have been rotated to keep the projection of the symmetry axis vertical in the image.
#     # This guarantees nothing about whether it tilts into the plane of the screen.
#     else:
#         if object2 == "8j1nA":
#             cmd.set_view(( \
#                 -0.931455255, -0.247361913, 0.266828001, \
#                 -0.333892405, 0.872502267, -0.356713563, \
#                 -0.144571185, -0.421355098, -0.895296514, \
#                 0.000000000, 0.000000000, -201.273925781, \
#                 116.961257935, 121.269165039, 142.889236450, \
#                 199.273925781, 243.273941040, -20.000000000))

#         elif object2 == "8hbwA" or object2 == "8g8wA":
#             cmd.set_view(( \
#                 -0.139716908, -0.225761563, 0.964110076, \
#                 -0.346182168, 0.923352599, 0.166048065, \
#                 -0.927702844, -0.310556084, -0.207162619, \
#                 0.000001798, -0.000004835, -205.644180298, \
#                 116.410858154, 121.141448975, 143.195556641, \
#                 202.644180298, 249.813613892, -20.000000000))

#         elif object2 == "8hbvA":
#             cmd.set_view(( \
#                 -0.313199431, -0.430462301, -0.846525192, \
#                 0.231343210, 0.829937279, -0.507619739, \
#                 0.921074748, -0.354825407, -0.160351112, \
#                 0.000000000, 0.000000000, -219.326431274, \
#                 116.949813843, 120.683715820, 141.931060791, \
#                 217.326431274, 261.326416016, -20.000000000))

#         cmd.clip("near", 40)

# #cmd.show("spheres", "object pseudo*")

# serial = -1
# if serial > 0:
#     #pdbdict = {"c": object1[0:4]}

#     cstring = ""
#     if cutaway == 1:
#         cstring="-cutaway"
#     elif cutaway == 2:
#         cstring="-surface"

#     cmd.png(f"/Users/jonathanborowsky/Documents/grabelab/aac1-ucp1/ucp1-{object2[0:4]}{cstring}-v{serial}", width=3000, height=3000, ray=1)


# if False:
#     # views for 8hbw

#     ### cut below here and paste into script ###
#     set_view (\
#          0.006601426,    0.633489966,    0.773721457,\
#         -0.187508881,    0.760797918,   -0.621309519,\
#         -0.982240438,   -0.140978307,    0.123807587,\
#          0.000000000,    0.000000000, -194.417327881,\
#        117.368011475,  120.961791992,  142.959869385,\
#        191.911804199,  235.522872925,  -20.000000000 )
#     ### cut above here and paste into script ###

#     ### cut below here and paste into script ###
#     set_view (\
#          0.067804530,    0.663157403,    0.745401084,\
#         -0.260148853,    0.733029783,   -0.628483236,\
#         -0.963184059,   -0.151301637,    0.222222015,\
#          0.000000000,    0.000000000, -194.417327881,\
#        117.368011475,  120.961791992,  142.959869385,\
#        192.129867554,  235.304809570,  -20.000000000 )
#     ### cut above here and paste into script ###

#     ### cut below here and paste into script ###
#     set_view (\
#          0.164676309,    0.644365788,    0.746774614,\
#         -0.316381186,    0.751614153,   -0.578773499,\
#         -0.934229434,   -0.140955016,    0.327638656,\
#          0.000000000,    0.000000000, -194.417327881,\
#        117.368011475,  120.961791992,  142.959869385,\
#        192.129867554,  235.304809570,  -20.000000000 )
#     ### cut above here and paste into script ###

#     ### cut below here and paste into script ###
#     set_view (\
#          0.437436014,    0.636432946,    0.635292888,\
#         -0.363335490,    0.771329463,   -0.522528052,\
#         -0.822577178,   -0.002253763,    0.568642557,\
#          0.000000000,    0.000000000, -219.326446533,\
#        117.368011475,  120.961791992,  142.959869385,\
#        216.714752197,  260.538116455,  -20.000000000 )
#     ### cut above here and paste into script ###

#     set_view( \
#         -0.976646483, -0.205126435, -0.063895635, \
#         -0.174888536, 0.931773245, -0.318134069, \
#         0.124791861, -0.299528629, -0.945890427, \
#         0.000019178, -0.000015303, -204.428924561, \
#         116.529212952, 120.836463928, 142.029006958, \
#         200.363311768, 248.894531250, -20.000000000)

#     #views for both 8j1n and 8hbw
#     ### cut below here and paste into script ###
#     set_view( \
#         -0.703500211, 0.454935998, 0.545996785, \
#         -0.084649146, 0.709156692, -0.699948668, \
#         -0.705633521, -0.538637161, -0.460383743, \
#         0.000138803, 0.000148659, -171.177749634, \
#         116.576019287, 119.829765320, 146.290237427, \
#         168.177749634, 200.977767944, -20.000000000)
#     ### cut above here and paste into script ###

#     set_view (\
#         -0.683030546,    0.485127181,    0.545996785,\
#         -0.053700220,    0.712169170,   -0.699948668,\
#         -0.728410840,   -0.507411242,   -0.460383743,\
#          0.000138803,    0.000148659, -199.770004272,\
#        116.576019287,  119.829765320,  146.290237427,\
#        196.770004272,  229.570022583,  -20.000000000 )