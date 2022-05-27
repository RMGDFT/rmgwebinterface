import streamlit as st
import py3Dmol
import numpy as np
from stmol import showmol
def view_xyz(crmg):
    xyz_data ="%d\n\n"%(len(crmg.atoms))
    for atom in  crmg.atoms:
        xyz_data += "%s  %f %f %f\n"%(atom[0], atom[1],atom[2],atom[3])


    xyz_view = py3Dmol.view()
    xyz_view.addModel(xyz_data, "xyz", 'mol')
                      
    scale = 0.18
    radius = 0.05
    xyz_view.setStyle({'sphere':{'colorscheme':'Jmol','scale':scale},
                       'stick':{'colorscheme':'Jmol', 'radius':radius}})
    xyz_view.zoomTo()
    xyz_view.render()
    showmol(xyz_view,height=500,width=800)
