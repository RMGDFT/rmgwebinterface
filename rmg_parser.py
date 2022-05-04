import streamlit as st
import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import CifFile
import subprocess
from utils import *
from uctools import *
from add_items import *

class rmg_interface():
    def vasp2cell(self, POSCAR=None):
        with open(POSCAR, "r") as f:
            all_lines = f.readlines()
        self.cell = CellData()
        self.cell.unit = "angstrom"
        self.cell.lengthscale = float(all_lines[1])   #alat
        self.ibrav = 0
        self.cell.latticevectors = []
        #  lines 3,4,5 are lattice vectors
        celldm= [] 
        for i in range(3):
            vec = all_lines[i+2].split()
            a0 = float(vec[0]) 
            a1 = float(vec[1]) 
            a2 = float(vec[2]) 
            self.cell.latticevectors.append([a0,a1,a2])
            celldm.append(sqrt(a0*a0+a1*a1+a2*a2))
        self.cell.a = celldm[0] * self.cell.lengthscale
        self.cell.b = celldm[1] * self.cell.lengthscale
        self.cell.c = celldm[2] * self.cell.lengthscale

        species_list = all_lines[5].split()
        num_atoms_spec =[int(a) for a in all_lines[6].split()]
        first_char = all_lines[7].strip()[0]
        self.atom_unit = ""
        atom_start_line = 8
        if first_char == "C" or first_char == "c":
            self.atom_unit = "Absolute"
        elif first_char == "D" or first_char == "d":
            self.atom_unit = "Cell Relative"
        else:
            first_char = all_lines[8].strip()[0]
            atom_start_line = 9
            if first_char == "C" or first_char == "c":
                self.atom_unit = "Absolute"
            elif first_char == "D" or first_char == "d":
                self.atom_unit = "Cell Relative"
            else:
                st.markdown("Error in POSCAR file")
        if self.atom_unit == "Absolute":
            scale = self.cell.lengthscale
        else:
            scale = 1.0
        self.atoms = []
        atom_line = atom_start_line -1    
        for isp in range(len(num_atoms_spec)):
            natom_per_sp = num_atoms_spec[isp]
            for n in range(natom_per_sp):
                atom_line += 1
                coor_str = all_lines[atom_line].split()
                self.atoms.append([species_list[isp], float(coor_str[0]) * scale, float(coor_str[1]) * scale, float(coor_str[2]) * scale]) 



    def xyz2cell(self, xyz_file=None): 
        #################################################################
        # Open and read xyz file
        self.cell = CellData()
        self.cell.unit = "angstrom"
        self.atom_unit = "Absolute"
        self.cell.lengthscale = 1.0

        with open(xyz_file, "r") as f:
            all_lines = f.readlines()
        num_atoms = int(all_lines[0])
        self.atoms = []
        for i in range(num_atoms):
            line = all_lines[i+2].split()
            x = float(line[1]) * self.cell.lengthscale
            y = float(line[2]) * self.cell.lengthscale
            z = float(line[3]) * self.cell.lengthscale
            self.atoms.append([line[0], x,y,z ] )
        ##########
        #  for xyz file, read-in the lattice information
        ###########
        # bounding box of the xyz atoms
        x_max = max(self.atoms, key=lambda x:x[1])[1]
        x_min = min(self.atoms, key=lambda x:x[1])[1]
        y_max = max(self.atoms, key=lambda x:x[2])[2]
        y_min = min(self.atoms, key=lambda x:x[2])[2]
        z_max = max(self.atoms, key=lambda x:x[3])[3]
        z_min = min(self.atoms, key=lambda x:x[3])[3]

        bounding_box = [x_min,x_max, y_min, y_max, z_min, z_max]
        (ibrav, a, b, c,latticevectors) = add_lattice(bounding_box)
        self.cell = CellData()
        self.ibrav = ibrav
        self.cell.a = a
        self.cell.b = b
        self.cell.c = c
        self.cell.latticevectors = latticevectors

    def cif2cell(self, cif_file=None): 
        #################################################################
        # Open and read CIF file
        cf = CifFile.ReadCif(cif_file)

        self.atom_unit = "Absolute"
        ##############################################
        # Get blocks
        cfkeys = list(cf.keys())
        cb = cf.get(cfkeys[0])
        # Get reference data
        ref = ReferenceData()
        ref.getFromCIF(cb)
        # Get cell data
        cd = CellData()
        cd.force = True
        cd.getFromCIF(cb)


        ##############################################
        # Generate cell
        if self.reducetoprim:
            cd.primitive()
        else:
            cd.conventional()

        inputcell = copy.copy(cd)

        self.cell = cd

        self.cell.newunit("bohr")

        self.ibrav = 0
        t = LatticeMatrix(self.cell.latticevectors)
        for i in range(3):
            for j in range(3):
                t[i][j] = self.cell.latticevectors[i][j]*self.cell.lengthscale
        ortho = abs(self.cell.a - t[0][0]) < 1.0e-5 
        ortho = ortho and abs(self.cell.b - t[1][1]) < 1.0e-5
        ortho = ortho and abs(self.cell.c - t[2][2]) < 1.0e-5

        if ortho: self.ibrav = 8
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        if system == 'cubic':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 1
                elif setting == 'F':
                    self.ibrav = 2
                elif setting == 'I':
                    self.ibrav = 3
            else:
                self.ibrav = 1
        if system == 'hexagonal':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 4
            #elif setting == 'R':
            #    return 5
        self.atoms = []
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(self.cell.latticevectors,b.position))
                for i in range(3):
                    t[i] = self.cell.lengthscale*t[i]
                self.atoms.append([b.spcstring(), t[0],t[1],t[2]])
        #print(self.atoms)

    def cell2rmg(self, mag):
        filestring = ""
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
        filestring += 'bravais_lattice_type="%s"  \n'%brav_type[self.ibrav]
        if self.cell.unit == "angstrom":
            filestring += 'crds_units = "Angstrom"  \n'
            filestring += 'lattice_units = "Angstrom"  \n'
        elif self.cell.unit == "bohr":
            filestring += 'crds_units = "Bohr"  \n'
            filestring += 'lattice_units = "Bohr"  \n'
        else:
            st.markdown("WARNING: unit = ",self.cell.unit)

        t = LatticeMatrix(self.cell.latticevectors)
        if self.ibrav !=0:
            filestring += 'a_length="%16.8f"  \n'%self.cell.a
            filestring += 'b_length="%16.8f"  \n'%self.cell.b
            filestring += 'c_length="%16.8f"  \n'%self.cell.c
        else:
            t = LatticeMatrix(self.cell.latticevectors)
            filestring += 'lattice_vector="  \n'
            for i in range(3):
                for j in range(3):
                    filestring += " %.12e "%(self.cell.latticevectors[i][j] * self.cell.lengthscale)
                filestring += '  \n'    
            filestring += '"  \n'

        filestring += 'atomic_coordinate_type = "%s"  \n'%self.atom_unit 
        filestring += 'atoms="  \n'
        atom_format = "%s  %.12e %.12e %.12e"
        iatom = 0
        for a in self.atoms:
            filestring += atom_format%(a[0],a[1], a[2], a[3])
            filestring += "  1 1 1 %6.2f %6.2f %6.2f  \n"%(mag[iatom][0],mag[iatom][1],mag[iatom][2])
            iatom += 1
        filestring += '"  \n'

        return filestring

    def __init__(self, filename, filetype):
        self.reducetoprim = True
        if not os.path.exists(filename):
            sys.stderr.write("***Error: The file "+filename+" could not be found.\n")
            sys.exit(2)
        if filetype == "cif":
            self.cif2cell(filename)
        elif filetype == "xyz":
            self.xyz2cell(filename)
        elif filetype == "vasp":
            self.vasp2cell(filename)

        # set up species list
        tmp = set([])
        for a in self.atoms:
                tmp.add(a[0])
        self.species = list(tmp)
        #self.rmginput = self.cell2rmg(mag)
