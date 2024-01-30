#!/usr/bin/python3
import numpy as np
import os
import sys
from rmg_parser import *
from add_items_text import *

rmg_branch_list = ["rmg base code", "rmg localized orbital module", "NEGF"]
rmg_branch = rmg_branch_list[0]
# p="rmg base code use delocalized real space grids as basis set, localized orbitals module optimized atom-centered localized orbitals as a basis set")

filetype_supported = ["cif", "xyz", "vasp", "QE"]
filetype =""


kdelt = 1000.0
cutoff = 100.0
pp_type_list = ["SG15(NC)", "GBRV-1.5(US)","Pseudo Dojo(NC)"]

if len(sys.argv) != 2:
    print("need to give a file including the initial files location and output directroy")
    exit()
with open(sys.argv[1], "r") as f:
   all_files = f.readlines()
for init_file in all_files:
   if(init_file.lstrip()[0] == "#"):
       continue

   line = init_file.split()
   filename = line[0]
   filetype = line[1]
   if not os.path.isdir(line[2]):
       os.mkdir(line[2])

   rmgfilename = line[2] +"/input"
   cutoff = float(line[3])
   kdelt = float(line[4])

   name_split = filename.split(".")
   
  # print("initial file name:", filename, "  filetype: ", filetype)
   
   crmg = rmg_interface(filename, filetype)
   
   description = filename
   rmginput_str = 'description="'+description+'"  \n'
   
   pp_type = pp_type_list[0]
   
   grid_lines = add_grid_text(crmg.cell, cutoff)
   pseudo_lines = add_pseudo_text(crmg.species, pp_type)
   kpoint_lines = add_kpoints_text(crmg.cell, kdelt)
   ctrl_lines = add_control_text()
   xc_lines = "\n"
   spin_lines, mag = add_spin_text(crmg.species_AFM, crmg.atoms, 0)
   
   rmginput_str += grid_lines
   rmginput_str += ctrl_lines
   rmginput_str += kpoint_lines
   rmginput_str += pseudo_lines
   rmginput_str += xc_lines
   
   if rmg_branch == "rmg base code":
       rmginput_str += crmg.cell2rmg(mag)
   else:
       print("only works for rmg base code now")
   
   
   with open(rmgfilename, "w") as f:
       f.write(rmginput_str)
