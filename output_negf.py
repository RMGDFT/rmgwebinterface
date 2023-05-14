import streamlit as st
import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import shutil
from zipfile import ZipFile
import tarfile

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

def output_negf(rmginput_str, crmg, a_lead1, a_lead2, a_center, nx_lead1, nx_lead2, nx_center, ny, nz, eq_lead, orbital_dict):
    
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
#drho_precond = "false"
#charge_pulay_Gspace = "false"
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

    for atom in crmg.atoms:
        if atom[1] < a_lead1:
            atoms_lead1.append(atom)
        elif atom[1] < a_lead1 + a_center:
            x = atom[1] - a_lead1
            atoms_center.append([atom[0], x, atom[2], atom[3]])
        else:
            x = atom[1] - a_lead1 - a_center
            atoms_lead2.append([atom[0], x, atom[2], atom[3]])


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
    if not os.path.isdir("NEGF_INPUTS/lead1"):
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
    input_center = atom_orbital_out(rmginput_str, a_center, nx_lead2, ny, nz, atoms_center, orbital_dict)
    with open(os.path.join("NEGF_INPUTS/center", "input"), "w") as f:
        f.write(input_center)

    a_3lead1 = 3.0 * a_lead1
    nx_3lead1 = 3*nx_lead1
    input_3lead1 = """
tart_mode_NEGF="111"
metalic="true"
num_blocks="3"
blocks_dim="%d %d %d"
potential_compass = "0 48 96 0 72 0 72"
chargedensity_compass = "0 48 96 0 72 0 72"
"""%(num_orb_lead1, num_orb_lead1, num_orb_lead1)

    input_3lead1 += atom_orbital_out(rmginput_str, a_3lead1, nx_3lead1, ny, nz, atoms_3lead1, orbital_dict)

    with open(os.path.join("NEGF_INPUTS/3lead_lead1", "input"), "w") as f:
        f.write(input_3lead1)

    a_3lead2 = 3.0 * a_lead2
    nx_3lead2 = 3*nx_lead2
    input_3lead2 = """
tart_mode_NEGF="111"
metalic="true"
num_blocks="3"
blocks_dim="%d %d %d"
potential_compass = "0 48 96 0 72 0 72"
chargedensity_compass = "0 48 96 0 72 0 72"
"""%(num_orb_lead2, num_orb_lead2, num_orb_lead2)

    input_3lead2 += atom_orbital_out(rmginput_str, a_3lead2, nx_3lead2, ny, nz, atoms_3lead2, orbital_dict)

    with open(os.path.join("NEGF_INPUTS/3lead_lead2", "input"), "w") as f:
        f.write(input_3lead2)


    with tarfile.open('NEGF_INPUTS.tar', "w:gz") as tar:
        tar.add('NEGF_INPUTS', arcname=os.path.basename('NEGF_INPUTS'))
    #shutil.make_archive('NEGF_INPUTS.zip', 'zip', root_dir='.', base_dir='NEGF_INPUTS')
    #with ZipFile('NEGF_INPUTS.zip', 'w') as zip_object:
    #    for folder_name, sub_folder, file_names in os.walk('NEGF_INPUTS'):
    #        for filename in file_names:
    #            file_path = os.path.join(folder_name, filename)
    #            zip_object.write(file_path, os.path.basename(file_path))
                
    with open('NEGF_INPUTS.tar', 'rb') as fp:

      st.download_button(
        label="Downlowd negf input files",
        data=fp,
        file_name = 'NEGF_INPUTS.tar',
        mime='application/zip'
      )

