import numpy as np
import mdtraj as md

def count_lipids(path, upper_leaflet_ref_query, lipid_resname, lipid_atom):
        
    trj = md.load(path)
    lipid_residues = [r for r in trj.top.residues if r.name == lipid_resname]
    lipid_resseqs = [r.resSeq for r in lipid_residues]
    
    ref_ind = trj.top.select(upper_leaflet_ref_query)
    ref_z_init = trj.xyz[0,ref_ind,2]
    
    lipid_inds = trj.top.select(f"resname {lipid_resname} and name {lipid_atom}")
    lipid_z_init = trj.xyz[0,lipid_inds,2]
    
    # plt.hist(lipid_z_init, bins=40)
    # plt.plot(ref_z_init, color="red")
    
    z_mean = np.mean(lipid_z_init)
    #print(z_mean)

    upperleaflet_rsq = []
    lowerleaflet_rsq = []
    
    upperleaflet = 0
    lowerleaflet = 0
    
    for li, lzi in zip(lipid_resseqs, lipid_z_init):
        if (lzi - z_mean)*np.sign(ref_z_init-z_mean) > 0:
            upperleaflet += 1
            upperleaflet_rsq.append(li)
        else:
            lowerleaflet += 1
            lowerleaflet_rsq.append(li)

    return upperleaflet_rsq, lowerleaflet_rsq, upperleaflet, lowerleaflet
