import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from rmg_parser import *
from add_items import *
st.title('RMG input User Interface')
st.markdown("Welcome to use rmg package! This interface will help you generate input file to run with our open-source rmg package. You can use one of the examples to show the options, you can also upload your own atomic structure file and choose options. Currently cif, xyz, and vasp formats are supported for loading the atomic structure information")
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: left}<style>',
        unsafe_allow_html=True)

uploaded_file = st.file_uploader("Uploda a file with cif, xyz or vasp format")
col1, col2 = st.columns(2)
example_ =  col1.checkbox("use an example ", False)
cif_or_xyz = "None"
if example_:
    cif_or_xyz = col2.radio("choose cif or xyz", ["None", "cif", "xyz", "vasp"])

filetype_supported = ["cif", "xyz", "vasp"]
filetype =""
if uploaded_file:
    if not os.path.isdir("tempDir"):
        os.mkdir("tempDir")
    with open(os.path.join("tempDir", uploaded_file.name), "wb") as f:
        f.write(uploaded_file.getbuffer())
    filename = "tempDir/"+uploaded_file.name
    name_split = filename.split(".")
    filext = ""
    if len(name_split) >1: filext = name_split[len(name_split)-1]
    if filext in filetype_supported: 
        filetype = filext
    else:
        filetype = st.radio("choose your filetype", ["None", "cif", "xyz", "vasp", "more is coming"])
elif cif_or_xyz == "cif":
    filename = "cifs/FeAs.cif"  
    filetype = "cif"
elif  cif_or_xyz == "xyz":
    filename = "xyz_files/C60.xyz"  
    filetype = "xyz"
elif cif_or_xyz == "vasp":
    filename = "vasp_examples/POSCAR_BN_cart"
    filetype = "vasp"
else:
    st.markdown("upload a file or choose an example")

if not (filetype in filetype_supported):
    st.markdown(filetype +  " filetype not programed")
else:
  st.markdown("File used: "+filename)
  description = st.text_input("description", value="description of the input file")
  rmginput_str = 'description="'+description+'"  \n'

  st.subheader("BASIC OPTIONS")
  crmg = rmg_interface(filename, filetype)
  grid_lines = add_grid(crmg.cell)
  pseudo_lines = add_pseudo(crmg.species)
  kpoint_lines = add_kpoints(crmg.cell)
  spin_lines, mag = add_spin(crmg.species, crmg.atoms)
  ctrl_lines = add_control()
  st.subheader("COMMONLY USED OPTIONS")
  scf_lines = add_scf()
  mixing_lines = add_mixing()
  xc_lines = add_xc(crmg.species)
  qmcpack_lines = add_qmcpack()
  IO_lines = add_IOctrl()

  st.subheader("PERFORMANCE TUNING OPTIONS")
  misc_lines = add_misc()

      
  rmginput_str += grid_lines
  rmginput_str += scf_lines
  rmginput_str += mixing_lines
  rmginput_str += xc_lines
  rmginput_str += qmcpack_lines
  rmginput_str += ctrl_lines
  rmginput_str += kpoint_lines
  rmginput_str += pseudo_lines
  rmginput_str += IO_lines
  rmginput_str += spin_lines
  rmginput_str += misc_lines

  rmginput_str += crmg.cell2rmg(mag)
  rmgfilename = os.path.basename(filename).split(".")[0] +".rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data=rmginput_str,
     file_name = rmgfilename)
  show_rmginput = st.checkbox("show the generated rmg input file", False)
  if show_rmginput:
    st.markdown(rmginput_str)

