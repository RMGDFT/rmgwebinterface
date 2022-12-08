import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from rmg_parser import *
from add_items import *
from view_xyz import *
from make_supercell import *
#from folder_picker import *
st.title('RMG Input File Generator')
st.markdown("Welcome to the RMG Input File Generator! RMG is an Open Source computer code for electronic structure calculations and modeling of materials and molecules. It is based on density functional theory and uses real space basis and pseudopotentials. Designed for scalability it has been run successfully on systems with thousands of nodes and hundreds of thousands of CPU cores. It runs on Linux/UNIX, Windows and OS X. This interface will help you generate input files to use with the rmg package. You can use one of the provided examples to examine the input options or you can upload your own atomic structure file and select options. Currently cif, xyz, and vasp formats are supported for loading the atomic structure information")
st.markdown("Without any changes the default calculation will use Norm-conserving pseudopotentials to perform an electronic quench at the gamma point.")
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: left}<style>',
        unsafe_allow_html=True)
rmg_branch = st.radio("Generate input for ", ["rmg base code", "rmg localized orbital module"], help="rmg base code use delocalized real space grids as basis set, localized orbitals module optimized atom-centered localized orbitals as a basis set")

uploaded_file = st.file_uploader("Begin by uploading a file in CIF, XYZ, or VASP(v5+) POSCAR, or QuantumEspresso format or using an example.")
col1, col2 = st.columns(2)
example_ =  col1.checkbox("use an example ", False)
cif_or_xyz = "None"
if example_:
    cif_or_xyz = col2.radio("choose cif or xyz", ["None", "cif", "xyz", "vasp", "QE"])

filetype_supported = ["cif", "xyz", "vasp", "QE"]
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
    elif "POSCAR" in uploaded_file.name:
        filetype = "vasp"
    else:
        filetype = st.radio("choose your filetype", ["None", "cif", "xyz", "vasp", "QE", "more is coming"])
elif cif_or_xyz == "cif":
    filename = "cifs/FeAs.cif"  
    filetype = "cif"
elif  cif_or_xyz == "xyz":
    filename = "example_C60/C60.xyz"  
    filetype = "xyz"
elif cif_or_xyz == "vasp":
    filename = "vasp_examples/POSCAR_BN_cart"
    filetype = "vasp"
elif cif_or_xyz == "QE":
    filename = "QE_examples/diamond.scf.in"
    filetype = "QE"
else:
    st.markdown("upload a file or choose an example")

if not (filetype in filetype_supported):
    st.markdown(filetype +  " filetype not programed")
else:
  st.markdown("File used: "+filename)

  desc_text = ""
  if filetype == "vasp":
      with open(filename, "r") as f:
          desc_text = f.readline()
  elif uploaded_file:
      desc_text = uploaded_file.name
  else:
      desc_text = filename;

  crmg = rmg_interface(filename, filetype)
 
  #folder_picker() 
  supercell = st.checkbox("make a supercell", False)
  supercell_name = ""
  if supercell:
      cell_mult = st.text_input("supercell expansion", "1  1  1")
      make_supercell(crmg, cell_mult)
      supercell_name=cell_mult.replace(" ","")

  
  view_str = st.checkbox("view atomic structure", False)
  if view_str:
      view_xyz(crmg)

  description = st.text_input("description", value=desc_text)
  rmginput_str = 'description="'+description+'"  \n'

  st.subheader("BASIC OPTIONS")
  grid_lines = add_grid(crmg.cell)
  if rmg_branch == "rmg localized orbital module":
     orbital_dict = add_orbital_info(crmg.species)
  pseudo_lines = add_pseudo(crmg.species)
  kpoint_lines = add_kpoints(crmg.cell)
  ctrl_lines = add_control()
  xc_lines = add_xc(crmg.species)
  spin_lines, mag = add_spin(crmg.species_AFM, crmg.atoms)
  st.subheader("COMMONLY USED OPTIONS")
  scf_lines=""
  mixing_lines = ""
  qmcpack_lines = ""
  IO_lines=""
  opts_more = st.checkbox("check the box for more options", False)
  if opts_more:
      scf_lines = add_scf()
      mixing_lines = add_mixing()
      qmcpack_lines = add_qmcpack()
      IO_lines = add_IOctrl()

  misc_lines = ""

  st.subheader("ADVANCED OPTIONS")
  misc_more = st.checkbox("check the box for other options", False)
   
  if misc_more:
    misc_lines = add_misc()

      
  rmginput_str += grid_lines
  rmginput_str += ctrl_lines
  rmginput_str += kpoint_lines
  rmginput_str += pseudo_lines
  rmginput_str += xc_lines
  rmginput_str += scf_lines
  rmginput_str += mixing_lines
  rmginput_str += qmcpack_lines
  rmginput_str += IO_lines
  rmginput_str += spin_lines
  rmginput_str += misc_lines

    
  if rmg_branch == "rmg base code":
    rmginput_str += crmg.cell2rmg(mag)
    rmgfilename = os.path.basename(filename).split(".")[0] + supercell_name+".rmg"
  else:
      rmginput_str += crmg.cell2rmg_on(mag, orbital_dict)
      rmgfilename = os.path.basename(filename).split(".")[0] + supercell_name+"_on.rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data=rmginput_str,
     file_name = rmgfilename)
  show_rmginput = st.checkbox("show the generated rmg input file", False)
  if show_rmginput:
    st.markdown(rmginput_str)

