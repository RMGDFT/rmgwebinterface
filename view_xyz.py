import streamlit as st
import py3Dmol
import numpy as np
from stmol import showmol
def view_xyz(crmg):
    xyz_data ="%d\n\n"%(len(crmg.atoms))
    
    if crmg.atom_unit == "Cell Relative":
        for atom in crmg.atoms:
            x = atom[1] * crmg.cell.latticevectors[0][0]
            x += atom[2] * crmg.cell.latticevectors[1][0]
            x += atom[3] * crmg.cell.latticevectors[2][0]
            y = atom[1] * crmg.cell.latticevectors[0][1]
            y += atom[2] * crmg.cell.latticevectors[1][1]
            y += atom[3] * crmg.cell.latticevectors[2][1]
            z = atom[1] * crmg.cell.latticevectors[0][2]
            z += atom[2] * crmg.cell.latticevectors[1][2]
            z += atom[3] * crmg.cell.latticevectors[2][2]
            x *= crmg.cell.lengthscale
            y *= crmg.cell.lengthscale
            z *= crmg.cell.lengthscale

            xyz_data += "%s  %f %f %f\n"%(atom[0], x,y,z)
    else:
        for atom in  crmg.atoms:
            xyz_data += "%s  %f %f %f\n"%(atom[0], atom[1] ,atom[2],atom[3])

    xyzview = py3Dmol.view(query='pdb:1A2C')
    xyzview.setStyle({'cartoon':{'color':'spectrum'}})
    showmol(xyzview, height = 500,width=800)
#   xyz_view = py3Dmol.view()
#   xyz_view.addModel(xyz_data, "xyz", 'mol')
#                     
#   scale = 0.18
#   radius = 0.05
#   xyz_view.setStyle({'sphere':{'colorscheme':'Jmol','scale':scale},
#                      'stick':{'colorscheme':'Jmol', 'radius':radius}})
#   #xyz_view.zoomTo()
#   #xyz_view.render()
#   showmol(xyz_view,height=500,width=800)
