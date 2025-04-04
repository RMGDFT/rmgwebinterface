try:
    import streamlit as st
    from add_items import *
except:
    from add_items_text import *
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

def token_to_value(all_lines, token):
    for line in all_lines:
        if(token in line):
            t_pos = line.find(token)
            v_pos_b = line.find("=", t_pos)
            v_pos_e = line.find(",", t_pos)
            return line[v_pos_b+1:v_pos_e]
    return None    
class rmg_interface():
    def qe2cell(self, qefile=None):
        with open(qefile, "r") as f:
            all_lines = f.readlines()
        self.cell = CellData()
        ibrav_str = token_to_value(all_lines, "ibrav")
        if ibrav_str == None:
            print("no ibrav from QE file")
            
        self.ibrav = int(ibrav_str)
        self.cell.latticevectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        celldm= [-1,-1,-1,-1,-1,-1] 
        for i in range(6):
            token = "celldm("+str(i+1)+")"
            celldm_str = token_to_value(all_lines, token)
            if celldm_str != None:
                celldm[i] = float(celldm_str)
        if self.ibrav == 1:
            self.cell.a = celldm[0] 
            self.cell.b = celldm[0] 
            self.cell.c = celldm[0]
            self.cell.latticevectors[0][0] = self.cell.a
            self.cell.latticevectors[1][1] = self.cell.b
            self.cell.latticevectors[2][2] = self.cell.c
        elif self.ibrav == 2:
            self.cell.a = celldm[0] 
            self.cell.b = celldm[0] 
            self.cell.c = celldm[0]
            self.cell.latticevectors[0][0] = 0.5 * celldm[0]
            self.cell.latticevectors[0][1] = 0.5 * celldm[0]
            self.cell.latticevectors[0][2] = 0.0
            self.cell.latticevectors[1][0] = 0.0
            self.cell.latticevectors[1][1] = 0.5 * celldm[0]
            self.cell.latticevectors[1][2] = 0.5 * celldm[0]
            self.cell.latticevectors[2][0] = 0.5 * celldm[0]
            self.cell.latticevectors[2][1] = 0.0
            self.cell.latticevectors[2][2] = 0.5 * celldm[0]

        elif self.ibrav == 3:
            self.cell.a = celldm[0] 
            self.cell.b = celldm[0] 
            self.cell.c = celldm[0]
            self.cell.latticevectors[0][0] = 0.5 * celldm[0]
            self.cell.latticevectors[0][1] = 0.5 * celldm[0]
            self.cell.latticevectors[0][2] = -0.5 * celldm[0]
            self.cell.latticevectors[1][0] = -0.5 * celldm[0]
            self.cell.latticevectors[1][1] = 0.5 * celldm[0]
            self.cell.latticevectors[1][2] = 0.5 * celldm[0]
            self.cell.latticevectors[2][0] = 0.5 * celldm[0]
            self.cell.latticevectors[2][1] = -0.5 * celldm[0]
            self.cell.latticevectors[2][2] = 0.5 * celldm[0]
        elif self.ibrav == 4:
            self.cell.a = celldm[0] 
            self.cell.b = celldm[0] 
            self.cell.c = celldm[2] * celldm[0]


            self.cell.latticevectors[0][0] = celldm[0]
            self.cell.latticevectors[1][0] = -0.5*celldm[0]
            self.cell.latticevectors[1][1] = 0.5 * sqrt(3.0) * celldm[0]
            self.cell.latticevectors[2][2] = self.cell.c
        elif self.ibrav == 8:
            self.cell.a = celldm[0] 
            self.cell.b = celldm[1] * celldm[0]
            self.cell.c = celldm[2] * celldm[0]
            self.cell.latticevectors[0][0] = self.cell.a
            self.cell.latticevectors[1][1] = self.cell.b
            self.cell.latticevectors[2][2] = self.cell.c
        elif self.ibrav == 0:
            cell_line_index = 0
            while i < len(all_lines):
                if "CELL_PARAMETERS" in all_lines[i]:
                    cell_line_index = i
                    break
                i = i+1
            cell_para_unit = all_lines[cell_line_index].replace("{"," ").replace("}"," ").split()[1]
            scale = 1.0
            if cell_para_unit == "alat":
                scale = celldm[0]
            elif cell_para_unit == "angstrom":
                scale = 1.0/0.529177
            elif cell_para_unit == "bohr":
                scale = 1.0
            for i in range(3):
                ai = all_lines[cell_line_index + i + 1].split()
                self.cell.latticevectors[i] = [float(ai[0]) *scale,float(ai[1]) *scale,float(ai[2]) *scale]
            self.cell.a = 0.0
            self.cell.b = 0.0
            self.cell.c = 0.0
            for i in range(3):
                self.cell.a += self.cell.latticevectors[0][i] * self.cell.latticevectors[0][i] 
                self.cell.b += self.cell.latticevectors[1][i] * self.cell.latticevectors[1][i] 
                self.cell.c += self.cell.latticevectors[2][i] * self.cell.latticevectors[2][i] 
            self.cell.a = sqrt(self.cell.a)
            self.cell.b = sqrt(self.cell.b)
            self.cell.c = sqrt(self.cell.c)
        else:
            print("ibrav not programedi for ibrav = ", self.ibrav)


        self.cell.unit = "bohr"
        num_atoms = int(token_to_value(all_lines, "nat"))
        atom_line_index = 0
        while i < len(all_lines):
            if "ATOMIC_POSITIONS" in all_lines[i]:
                atom_line_index = i
                break
            i = i+1
        atom_unit = all_lines[atom_line_index].replace("{"," ").replace("}"," ").split()
        self.atom_unit = "Absolute"
        self.atoms = []
        scale = 1.0
        if(len(atom_unit) == 1):
            self.atom_unit = "Cell Relative"
        elif atom_unit[1] == "alat":
            scale = celldm[0]
        elif atom_unit[1] == "angstrom":
            scale = 1.0/0.529177
        elif atom_unit[1] == "bohr":
            scale = 1.0
        elif atom_unit[1] == "crystal":
            self.atom_unit = "Cell Relative"
        else:
            print("crystal_sg not programed")

        for i in range(num_atoms):
            oneatom = all_lines[atom_line_index + i + 1].split()
            self.atoms.append([oneatom[0], float(oneatom[1]) * scale, float(oneatom[2]) * scale, float(oneatom[3]) * scale ])

        


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
        ibrav = 0
        a = 1.0
        b = 1.0
        c = 1.0
        need_lattice = True
        latticevectors = []
        if "Lattice" in all_lines[1]:
           lat = all_lines[1].split('"')[1] 
           if len(lat.split()) == 9:
                a0 = float(lat.split()[0])
                a1 = float(lat.split()[1])
                a2 = float(lat.split()[2])
                latticevectors.append([a0, a1, a2])
                a = sqrt(a0*a0+a1*a1+a2*a2)
                a0 = float(lat.split()[3])
                a1 = float(lat.split()[4])
                a2 = float(lat.split()[5])
                latticevectors.append([a0, a1, a2])
                b = sqrt(a0*a0+a1*a1+a2*a2)
                a0 = float(lat.split()[6])
                a1 = float(lat.split()[7])
                a2 = float(lat.split()[8])
                latticevectors.append([a0, a1, a2])
                c = sqrt(a0*a0+a1*a1+a2*a2)

                need_lattice = False
                centered = False
        if need_lattice:
            # bounding box of the xyz atoms
            x_max = max(self.atoms, key=lambda x:x[1])[1]
            x_min = min(self.atoms, key=lambda x:x[1])[1]
            y_max = max(self.atoms, key=lambda x:x[2])[2]
            y_min = min(self.atoms, key=lambda x:x[2])[2]
            z_max = max(self.atoms, key=lambda x:x[3])[3]
            z_min = min(self.atoms, key=lambda x:x[3])[3]

            bounding_box = [x_min,x_max, y_min, y_max, z_min, z_max]
            (ibrav, a, b, c,latticevectors, centered) = add_lattice(bounding_box)
        self.cell = CellData()
        self.ibrav = ibrav
        self.cell.a = a
        self.cell.b = b
        self.cell.c = c
        self.cell.latticevectors = latticevectors
        if centered and ibrav == 8:
            xshift = (x_min + x_max -a) *0.5 
            yshift = (y_min + y_max -b) *0.5 
            zshift = (z_min + z_max -c) *0.5 
            for atom in self.atoms:
                atom[1] -= xshift
                atom[2] -= yshift
                atom[3] -= zshift


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
        filestring = "#****  LATTICE and ATOMS  ****   \n"
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
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]

            filestring += atom_format%(b,a[1], a[2], a[3])
            filestring += "  1 1 1 %6.2f %6.2f %6.2f  \n"%(mag[iatom][0],mag[iatom][1],mag[iatom][2])
            iatom += 1
        filestring += '"  \n'

        return filestring

    def cell2rmg_on(self, mag, orbital_dict):
        filestring = "#*** for ON and NEGF, alyways use multigrid level 2 *** \n"
        filestring += 'kohn_sham_mg_levels ="2" \n'
        filestring += '#*** non-local projector should be always localized for ON and NEGF***i \n'
        filestring += 'max_nlradius="6.0" \n\n'
        filestring += "#****  LATTICE and ATOMS  ****   \n"
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
        # force to use atomic unit 
        filestring += 'crds_units = "Bohr"  \n'
        filestring += 'lattice_units = "Bohr"  \n'

        if self.cell.unit == "angstrom":
            bohr= 1.88972612499359;
        elif self.cell.unit == "bohr":
            bohr = 1.0
        else:
            st.markdown("WARNING: unit = ",self.cell.unit)

        self.cell.lengthscale = self.cell.lengthscale * bohr
        t = LatticeMatrix(self.cell.latticevectors)
        if self.ibrav !=0:
            filestring += 'a_length="%16.8f"  \n'%(self.cell.a * bohr)
            filestring += 'b_length="%16.8f"  \n'%(self.cell.b * bohr)
            filestring += 'c_length="%16.8f"  \n'%(self.cell.c * bohr)
        else:
            t = LatticeMatrix(self.cell.latticevectors)
            filestring += 'lattice_vector="  \n'
            for i in range(3):
                for j in range(3):
                    filestring += " %.12e "%(self.cell.latticevectors[i][j] * self.cell.lengthscale)
                filestring += '  \n'    
            filestring += '"  \n'


        num_orb_tot = 0
        for i in range(len(self.atoms)):
            b = self.atoms[i][0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]
            num_orb_tot += orbital_dict[b][0]
            x = self.atoms[i][1] * bohr
            y = self.atoms[i][2] * bohr
            z = self.atoms[i][3] * bohr
            if self.atom_unit == "Cell Relative":
                x = self.atoms[i][1] * self.cell.latticevectors[0][0]
                x += self.atoms[i][2] * self.cell.latticevectors[1][0]
                x += self.atoms[i][3] * self.cell.latticevectors[2][0]
                x *= self.cell.lengthscale
                y = self.atoms[i][1] * self.cell.latticevectors[0][1]
                y += self.atoms[i][2] * self.cell.latticevectors[1][1]
                y += self.atoms[i][3] * self.cell.latticevectors[2][1]
                y *= self.cell.lengthscale
                z = self.atoms[i][1] * self.cell.latticevectors[0][2]
                z += self.atoms[i][2] * self.cell.latticevectors[1][2]
                z += self.atoms[i][3] * self.cell.latticevectors[2][2]
                z *= self.cell.lengthscale

            self.atoms[i] = [self.atoms[i][0], x,y,z] + mag[i]
        self.atoms.sort(key=lambda x:x[1])

        filestring += 'number_of_orbitals = "%d"  \n'%num_orb_tot
        filestring += 'number_of_atoms = "%d"  \n'%(len(self.atoms))
        filestring += 'number_of_species = "%d"  \n'%(len(self.species))

        filestring += 'atomic_coordinate_type = "Absolute"  \n'
        atom_format = "%s  %.12e %.12e %.12e"
        orbital_format = "%d  %.12e %.12e %.12e %7.4f   1  1  \n"
        filestring += 'atoms="  \n'
        for a in self.atoms:
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]
            filestring += atom_format%(b,a[1], a[2], a[3])
            filestring += "  1 1 1 %6.2f %6.2f %6.2f  \n"%(a[4], a[5], a[6])
        filestring += '"  \n'

        filestring += '#Orbitalsi: number per center, center coordinates and radious  \n'
        filestring += 'orbitals = "  \n'
        for a in self.atoms:
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]
            filestring += orbital_format%(orbital_dict[b][0], a[1], a[2],a[3],orbital_dict[b][1])

        filestring += '"  \n'
        filestring += 'LocalizedOrbitalLayout = "Projection"  \n'


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
        elif filetype == "QE":
            self.qe2cell(filename)

        # set up species list
        tmp = set([])
        tmp_AFM = set([])
        for a in self.atoms:
            b = a[0]
            if b[len(b) -1].isdigit():
                b = b[:len(b)-1]
            tmp_AFM.add(a[0])
            tmp.add(b)
        self.species = list(tmp)
        self.species_AFM = list(tmp_AFM)
        #self.rmginput = self.cell2rmg(mag)
