import count_lipids

# CELL 3
# GENERATE MDTRAJ QUERY FOR SPECIFIED (RESIDUE, LEAFLET, PROTEIN) COMBINATION
##########################################################################################

def get_ligand_query(protein, residue, leaflet, refpath):
    #parameters
    names = {"ARAN":["C1", "C5", "C10", "C15", "C20"], "POPC":["P"]}
    namequery = " or ".join(["name "+ name for name in names[residue]])

    #select lipids or fatty acids
    if residue == "POPC":
        init_rsqs = {"aac1":1, "ucp1":9}
        upper_leaflet_ref_query = f"resSeq {init_rsqs[protein]} and name CA"
        upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet = count_lipids.count_lipids(refpath, upper_leaflet_ref_query, residue, names[residue][0])
        if leaflet == "upper":
            leaflet_rsq = upperleaflet_rsq
            nleaflet = upperleaflet
        elif leaflet == "lower":
            leaflet_rsq = lowerleaflet_rsq
            nleaflet = lowerleaflet
        else: 
            raise ValueError("leaflet must be upper or lower")
        
        c_l_query = " or ".join(["resSeq "+str(i) for i in leaflet_rsq])
        ligand_query = f"resname {residue} and {namequery} and ({c_l_query})"

    elif residue == "ARAN":
        c_aran = {"aac1":[309,310,311,312,313,314,315,316], "ucp1":[527,528,529,530,531,532,533,534]}
        invstring = ""
        if leaflet == "lower":
            invstring = " not "
        c_aran_query = " or ".join(["resSeq "+str(i) for i in c_aran[protein]])
        ligand_query = f"resname {residue} and ({namequery}) and {invstring}({c_aran_query})"
        nleaflet = 8 #number of ARAN in each leaflet

    else:
        raise ValueError("residue must be POPC or ARAN")

    return names[residue], f"({ligand_query})", nleaflet