import streamlit as st
import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import tarfile
from config_part import *



def LCR_file_output(lead1, center, lead2):

    if center.dirname == lead1.dirname:
        dir_3lead1 = lead1.dirname
        dir_3lead2 = lead2.dirname
    else:
        dir_3lead1 = "3lead_lead1"
        dir_3lead2 = "3lead_lead2"

    lcr1 = """
        #  orbitals are read from
../%s/Waves/wave.out

#  charge density and potentials are read from
../%s/Waves/wave.out
"""% (lead1.dirname, dir_3lead1)

    lcr1 += """
#   lcr[].NX_GRID, NY_GRID, NZ_GRID
%d  %d  %d

#   starting grid point in the orignal grid
0  0  0

#   ending grid point in the orignal grid
%d  %d  %d

#   starting grid point in the NEGF globla grid
0  0  0
"""%(lead1.nx, lead1.ny, lead1.nz, lead1.nx, lead1.ny, lead1.nz )

    lcr1 += """
# lcr[].num_ions  number of ions in a conductor or a lead
%d

#lcr[].state_begin, state_middle, state_end state_middle is dummy
0  0  %d
"""%(lead1.num_atoms, lead1.num_orb)

    lcr1 += """
#lcr[].ion_begin
  0

#lcr[].bias (eV)  dummy
0.0       0.0       0.0

# a_length of this supercell (same as ON2 cal),
# the starting x-coordinate (in the global system)
%f  0.0

# b_length of this supercell (same as ON2 cal),
# the starting y-coordinate (in the global system)
%f  0.0
"""%(lead1.a, lead1.b)

    lcr2 = """
        #  orbitals are read from
../%s/Waves/wave.out

#  charge density and potentials are read from
../%s/Waves/wave.out
"""% (lead2.dirname, dir_3lead2)

    lcr2 += """
#   lcr[].NX_GRID, NY_GRID, NZ_GRID
%d  %d  %d

#   starting grid point in the orignal grid
0  0  0

#   ending grid point in the orignal grid
%d  %d  %d

#   starting grid point in the NEGF globla grid
0  0  0
"""%(lead2.nx, lead2.ny, lead2.nz, lead2.nx, lead2.ny, lead2.nz )

    lcr2 += """
# lcr[].num_ions  number of ions in a conductor or a lead
%d

#lcr[].state_begin, state_middle, state_end state_middle is dummy
0  0  %d
"""%(lead2.num_atoms, lead2.num_orb)

    lcr2 += """
#lcr[].ion_begin
  0

#lcr[].bias (eV)  dummy
0.0       0.0       0.0

# a_length of this supercell (same as ON2 cal),
# the starting x-coordinate (in the global system)
%f  %f

# b_length of this supercell (same as ON2 cal),
# the starting y-coordinate (in the global system)
%f  0.0
"""%(lead2.a, lead1.a + center.a, lead2.b)

    lcr0 = """
        #  orbitals are read from
../%s/Waves/wave.out

#  charge density and potentials are read from
../%s/Waves/wave.out
"""% (center.dirname, center.dirname)

    lcr0 += """
#   lcr[].NX_GRID, NY_GRID, NZ_GRID
%d  %d  %d

#   starting grid point in the orignal grid
0  0  0

#   ending grid point in the orignal grid
%d  %d  %d

#   starting grid point in the NEGF globla grid
0  0  0
"""%(center.nx, center.ny, center.nz, center.nx, center.ny, center.nz )

    lcr0 += """
# lcr[].num_ions  number of ions in a conductor or a lead
%d

#lcr[].state_begin, state_middle, state_end state_middle is dummy
0  0  %d
"""%(center.num_atoms, center.num_orb)

    lcr0 += """
#lcr[].ion_begin
  0

#lcr[].bias (eV)  dummy
0.0       0.0       0.0

# a_length of this supercell (same as ON2 cal),
# the starting x-coordinate (in the global system)
%f  %f

# b_length of this supercell (same as ON2 cal),
# the starting y-coordinate (in the global system)
%f  0.0
"""%(center.a, lead1.a, center.b)
    trans = """
#num of probes and which block the probe attached to
2  0  2

#num of sub systems and their atomic order in NEGF calculations
#need to edit for multi-probe calclations
3  1  0  2
# num_probe_potential_window & their ranges[in the order of lead 1, 2, 3 ...]
#they are dummy for 2-probe calculations and 
#need to edit for multi-probe calclations
4  40 120  40 120  40 120  40 120

# num_dos_axis_window for integration along x,y,z axis
# (required for mode 200)
3  0 110  0 72 0 72

#ncircle 
 50 


#nmax_gq1 
 50 


#nmax_gq2 
 128 


#Enery Low Bound 
 -30.0 


#KT 
 0.025 


#GAMMA 
 0.5 


#DELTA2 
 2.05 


#DELTA 
 1.0E-6 


#Charge density Pulay order 
 5 


#Charge density Pulay refresh steps 
 100 


#Charge density Pulay beta 
 0.5 


#processor grids for block matrix operation  
1 
1
"""

    cond_input = """
#  input for conductance calculations 

%d %d %d  nC nL nR 
3 %d %d %d  num_blocks, block_dim 
-2.0 2.0 401 Emin,Emax, E_points 
0.005  small imaginary part 
0.05  kbt 
2001 1 1   
1   number of conductance curve from lead i to lead j
1  2  lead i to lead j
"""%(lead1.num_orb+lead2.num_orb+center.nx, lead1.num_orb, lead2.num_orb, lead1.num_orb, center.num_orb, lead2.num_orb)
    return lcr0, lcr1, lcr2, trans, cond_input

def atom_orbital_out(rmginput_str, a, nx, ny, nz, atoms, orbital_dict):
    
    num_orb_tot = 0
    species_list = []
    for atom in atoms:
        sp = atom[0]
        num_orb_tot += orbital_dict[sp][0]
        species_list.append(sp)
    species = set(species_list)

    filestring = ""
    filestring += 'a_length="%16.8f"  \n'%(a)
    filestring += 'wavefunction_grid="%d %d %d"  \n'%(nx,ny,nz)
    filestring += 'number_of_orbitals = "%d"  \n'%num_orb_tot
    filestring += 'number_of_atoms = "%d"  \n'%(len(atoms))
    filestring += 'number_of_species = "%d"  \n'%(len(species))

    filestring += 'atomic_coordinate_type = "Absolute"  \n'
    atom_format = "%s  %.12e %.12e %.12e  1  \n"
    orbital_format = "%d  %.12e %.12e %.12e %7.4f   1  1  \n"
    filestring += 'atoms="  \n'
    for a in atoms:
        b = a[0]
        if b[len(b) -1].isdigit():
            b = b[:len(b)-1]
        filestring += atom_format%(b,a[1], a[2], a[3])
    filestring += '"  \n'

    filestring += '#Orbitalsi: number per center, center coordinates and radious  \n'
    filestring += 'orbitals = "  \n'
    for a in atoms:
        b = a[0]
        if b[len(b) -1].isdigit():
            b = b[:len(b)-1]
        filestring += orbital_format%(orbital_dict[b][0], a[1], a[2],a[3],orbital_dict[b][1])

    filestring += '"  \n'


    return rmginput_str + filestring

def output_negf(rmginput_str, crmg, a_lead1, a_lead2, a_center, nx_lead1, nx_lead2, nx_center, ny, nz, eq_lead, num_atoms_lead1, num_atoms_lead2, orbital_dict):
    
    if not eq_lead:
       st.markdown("WARNING: not programmed yet")
       return


    rmginput_str += "#****  LATTICE and ATOMS  ****   \n"
    #
    # some default input options
    brav_type = {
        0:"None",
        1:"Cubic Primitive",
        2:"Cubic Face Centered",
        3:"Cubic Body Centered",
        4:"Hexagonal Primitive",
        8:"Orthorhombic Primitive"
    }
    rmginput_str += 'bravais_lattice_type="%s"  \n'%brav_type[8]
    # force to use atomic unit 
    # **** SCF controls ****

    some_default ="""
max_scf_steps = "100"
freeze_orbital_step="80"
rms_convergence_criterion = "1e-7"

# **** Mixing controls ****
orbital_mixing = "0.1"
orbital_pulay_refresh="100"
sqrt_interpolation = "true"

charge_density_mixing = "0.1"
charge_mixing_type = "Pulay"
charge_pulay_order = "5"
charge_pulay_scale = "0.1"
charge_pulay_refresh = "20"
drho_precond = "false"
charge_pulay_Gspace = "false"
occupation_electron_temperature_eV = "0.1"
kohn_sham_mg_levels = "2"
"""


    rmginput_str += some_default
    rmginput_str += 'max_nlradius = "6.00000000"  \n'
    rmginput_str += 'crds_units = "Bohr"  \n'
    rmginput_str += 'lattice_units = "Bohr"  \n'
    rmginput_str += 'LocalizedOrbitalLayout = "Projection"  \n'
    rmginput_str += 'localize_localpp = "true"  \n'
    rmginput_str += 'localize_projectors = "true"  \n'


    if crmg.cell.unit == "angstrom":
        bohr= 1.88972612499359;
    elif crmg.cell.unit == "bohr":
        bohr = 1.0
    else:
        st.markdown("WARNING: unit = ",crmg.cell.unit)

    crmg.cell.lengthscale = crmg.cell.lengthscale * bohr
    rmginput_str += 'b_length="%16.8f"  \n'%(crmg.cell.b * bohr)
    rmginput_str += 'c_length="%16.8f"  \n'%(crmg.cell.c * bohr)



    num_orb_tot = 0
    for i in range(len(crmg.atoms)):
        b = crmg.atoms[i][0]
        if b[len(b) -1].isdigit():
            b = b[:len(b)-1]
        num_orb_tot += orbital_dict[b][0]
        x = crmg.atoms[i][1] * bohr
        y = crmg.atoms[i][2] * bohr
        z = crmg.atoms[i][3] * bohr
        if crmg.atom_unit == "Cell Relative":
            x = crmg.atoms[i][1] * crmg.cell.latticevectors[0][0]
            x += crmg.atoms[i][2] * crmg.cell.latticevectors[1][0]
            x += crmg.atoms[i][3] * crmg.cell.latticevectors[2][0]
            x *= crmg.cell.lengthscale
            y = crmg.atoms[i][1] * crmg.cell.latticevectors[0][1]
            y += crmg.atoms[i][2] * crmg.cell.latticevectors[1][1]
            y += crmg.atoms[i][3] * crmg.cell.latticevectors[2][1]
            y *= crmg.cell.lengthscale
            z = crmg.atoms[i][1] * crmg.cell.latticevectors[0][2]
            z += crmg.atoms[i][2] * crmg.cell.latticevectors[1][2]
            z += crmg.atoms[i][3] * crmg.cell.latticevectors[2][2]
            z *= crmg.cell.lengthscale

        crmg.atoms[i] = [crmg.atoms[i][0], x,y,z]
    crmg.atoms.sort(key=lambda x:x[1])

    a_lead1 *= bohr
    a_lead2 *= bohr
    a_center *= bohr
    atoms_lead1 = []
    atoms_lead2 = []
    atoms_center = []
    atoms_3lead1 = []
    atoms_3lead2 = []

    if crmg.atoms[0][1] < 0.0:
        for i in range(len(crmg.atoms)):
            crmg.atoms[i][1] += crmg.atoms[0][1]
            
    num_atoms_center = len(crmg.atoms) - num_atoms_lead1 - num_atoms_lead2

    atoms_lead1 = crmg.atoms[0:num_atoms_lead1]
    atoms_center = crmg.atoms[num_atoms_lead1:num_atoms_lead1 + num_atoms_center]
    atoms_lead2 = crmg.atoms[num_atoms_lead1+num_atoms_center:]

    for i in range(len(atoms_center)):
        atoms_center[i][1] -= a_lead1
    for i in range(len(atoms_lead2)):
        atoms_lead2[i][1] -= a_lead1 + a_center

    for i in range(3):
        for atom in atoms_lead1:
            x = atom[1] + i * a_lead1
            atoms_3lead1.append([atom[0], x, atom[2], atom[3]])
    for i in range(3):
        for atom in atoms_lead2:
            x = atom[1] + i * a_lead2
            atoms_3lead2.append([atom[0], x, atom[2], atom[3]])

    num_orb_lead1 = 0
    for atom in atoms_lead1:
        sp = atom[0]
        num_orb_lead1 += orbital_dict[sp][0]

    num_orb_lead2 = 0
    for atom in atoms_lead2:
        sp = atom[0]
        num_orb_lead2 += orbital_dict[sp][0]

    num_orb_center = 0
    for atom in atoms_center:
        sp = atom[0]
        num_orb_center += orbital_dict[sp][0]


    lead1 = config_part("lead1", a_lead1,crmg.cell.b * bohr, crmg.cell.c*bohr, nx_lead1, ny, nz, num_atoms_lead1, num_orb_lead1)
    lead2 = config_part("lead2", a_lead2,crmg.cell.b * bohr, crmg.cell.c*bohr, nx_lead2, ny, nz, num_atoms_lead2, num_orb_lead2)
    center = config_part("center", a_center,crmg.cell.b * bohr, crmg.cell.c*bohr, nx_center, ny, nz, num_atoms_center, num_orb_center)
    input_lead1 = rmginput_str
    input_lead2 = rmginput_str
    input_3lead1 = rmginput_str
    input_3lead2 = rmginput_str
    input_center = rmginput_str
    input_bias = rmginput_str


    if not os.path.isdir("NEGF_INPUTS"):
        os.mkdir("NEGF_INPUTS")
    if not os.path.isdir("NEGF_INPUTS/lead1"):
        os.mkdir("NEGF_INPUTS/lead1")
    if not os.path.isdir("NEGF_INPUTS/lead2"):
        os.mkdir("NEGF_INPUTS/lead2")
    if not os.path.isdir("NEGF_INPUTS/center"):
        os.mkdir("NEGF_INPUTS/center")
    if not os.path.isdir("NEGF_INPUTS/3lead_lead1"):
        os.mkdir("NEGF_INPUTS/3lead_lead1")
    if not os.path.isdir("NEGF_INPUTS/3lead_lead2"):
        os.mkdir("NEGF_INPUTS/3lead_lead2")
    if not os.path.isdir("NEGF_INPUTS/bias_0.0"):
        os.mkdir("NEGF_INPUTS/bias_0.0")

    input_lead1 = atom_orbital_out(rmginput_str, a_lead1, nx_lead1, ny, nz, atoms_lead1, orbital_dict)
    with open(os.path.join("NEGF_INPUTS/lead1", "input"), "w") as f:
        f.write(input_lead1)
    input_lead2 = atom_orbital_out(rmginput_str, a_lead2, nx_lead2, ny, nz, atoms_lead2, orbital_dict)
    with open(os.path.join("NEGF_INPUTS/lead2", "input"), "w") as f:
        f.write(input_lead2)
    input_center = atom_orbital_out(rmginput_str, a_center, nx_center, ny, nz, atoms_center, orbital_dict)
    with open(os.path.join("NEGF_INPUTS/center", "input"), "w") as f:
        f.write(input_center)

    a_3lead1 = 3.0 * a_lead1
    nx_3lead1 = 3*nx_lead1
    input_3lead1 = """
start_mode_NEGF="111"
metalic="true"
num_blocks="3"
blocks_dim="%d %d %d"
potential_compass = "0 48 96 0 72 0 72"
chargedensity_compass = "0 48 96 0 72 0 72"
"""%(num_orb_lead1, num_orb_lead1, num_orb_lead1)

    input_3lead1 += atom_orbital_out(rmginput_str, a_3lead1, nx_3lead1, ny, nz, atoms_3lead1, orbital_dict)

    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "input"), "w") as f:
        f.write(input_3lead1)

    lcr0, lcr1, lcr2, tran, cond = LCR_file_output(lead1, lead1, lead1)

    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "LCR.dat0"), "w") as f:
        f.write(lcr0)
    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "LCR.dat1"), "w") as f:
        f.write(lcr1)
    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "LCR.dat2"), "w") as f:
        f.write(lcr2)
    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "trans.in"), "w") as f:
        f.write(tran)
    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "cond.input"), "w") as f:
        f.write(cond)


    lcr0, lcr1, lcr2, tran, cond = LCR_file_output(lead2, lead2, lead2)

    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "LCR.dat0"), "w") as f:
        f.write(lcr0)
    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "LCR.dat1"), "w") as f:
        f.write(lcr1)
    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "LCR.dat2"), "w") as f:
        f.write(lcr2)
    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "trans.in"), "w") as f:
        f.write(tran)
    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "cond.input"), "w") as f:
        f.write(cond)

    lcr0, lcr1, lcr2, tran, cond = LCR_file_output(lead1, center, lead2)

    with open(os.path.join("NEGF_INPUTS/bias_0.0", "LCR.dat0"), "w") as f:
        f.write(lcr0)
    with open(os.path.join("NEGF_INPUTS/bias_0.0", "LCR.dat1"), "w") as f:
        f.write(lcr1)
    with open(os.path.join("NEGF_INPUTS/bias_0.0", "LCR.dat2"), "w") as f:
        f.write(lcr2)
    with open(os.path.join("NEGF_INPUTS/bias_0.0", "trans.in"), "w") as f:
        f.write(tran)
    with open(os.path.join("NEGF_INPUTS/bias_0.0", "cond.input"), "w") as f:
        f.write(cond)

    a_3lead2 = 3.0 * a_lead2
    nx_3lead2 = 3*nx_lead2
    input_3lead2 = """
start_mode_NEGF="111"
metalic="true"
num_blocks="3"
blocks_dim="%d %d %d"
potential_compass = "0 48 96 0 72 0 72"
chargedensity_compass = "0 48 96 0 72 0 72"
"""%(num_orb_lead2, num_orb_lead2, num_orb_lead2)

    input_3lead2 += atom_orbital_out(rmginput_str, a_3lead2, nx_3lead2, ny, nz, atoms_3lead2, orbital_dict)

    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "input"), "w") as f:
        f.write(input_3lead2)

    input_bias = """
start_mode_NEGF="112"
num_blocks="3"
blocks_dim="%d %d %d"
"""%(num_orb_lead1, num_orb_center, num_orb_lead2)

    input_bias+= """
potential_compass = "1 %d %d 0 %d 0 %d"
"""%(nx_lead1, nx_lead1 + nx_center, ny, nz)

    input_bias+= """
chargedensity_compass = "1 %d %d 0 %d 0 %d"
"""%(nx_lead1, nx_lead1 + nx_center, ny, nz)

    input_bias += atom_orbital_out(rmginput_str, a_lead1+a_center+a_lead2, nx_lead1 + nx_lead2 + nx_center, ny, nz, crmg.atoms, orbital_dict)

    with open(os.path.join("NEGF_INPUTS/bias_0.0", "input"), "w") as f:
        f.write(input_bias)

    input_bias = input_bias.replace('start_mode_NEGF="112"', 'start_mode_NEGF="110"')
    with open(os.path.join("NEGF_INPUTS/bias_0.0", "input.110"), "w") as f:
        f.write(input_bias)

    rmgfilename = 'NEGF_INPUTS'
    rmgfilename = st.text_input("output file name", rmgfilename)
    rmgfilename += '.tar'
    with tarfile.open(rmgfilename, "w:gz") as tar:
        tar.add('NEGF_INPUTS', arcname=os.path.basename('NEGF_INPUTS'))
    #shutil.make_archive('NEGF_INPUTS.zip', 'zip', root_dir='.', base_dir='NEGF_INPUTS')
    #with ZipFile('NEGF_INPUTS.zip', 'w') as zip_object:
    #    for folder_name, sub_folder, file_names in os.walk('NEGF_INPUTS'):
    #        for filename in file_names:
    #            file_path = os.path.join(folder_name, filename)
    #            zip_object.write(file_path, os.path.basename(file_path))
                
    with open(rmgfilename, 'rb') as fp:

      st.download_button(
        label="Downlowd negf input files",
        data=fp,
        file_name = rmgfilename,
        mime='application/zip'
      )

