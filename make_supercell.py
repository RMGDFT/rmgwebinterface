import streamlit as st
import numpy as np
def make_supercell(crmg, cell_mult):
    cells = [int(i) for i in cell_mult.split()]
    crmg.ibrav = 0
    crmg.cell.a = crmg.cell.a * cells[0]
    crmg.cell.b = crmg.cell.b * cells[1]
    crmg.cell.c = crmg.cell.c * cells[2]
    atoms_super = []

    if crmg.atom_unit == "Absolute":
       for i in range(cells[0]):
        for j in range(cells[1]):
            for k in range(cells[2]):
                x = i* crmg.cell.latticevectors[0][0] + j* crmg.cell.latticevectors[1][0] + k* crmg.cell.latticevectors[2][0] 
                y =   i* crmg.cell.latticevectors[0][1] + j* crmg.cell.latticevectors[1][1] + k* crmg.cell.latticevectors[2][1] 
                z =   i* crmg.cell.latticevectors[0][2] + j* crmg.cell.latticevectors[1][2] + k* crmg.cell.latticevectors[2][2] 

                x = x * crmg.cell.lengthscale
                y = y * crmg.cell.lengthscale
                z = z * crmg.cell.lengthscale
                for atom in  crmg.atoms:
                    new_atom = [atom[0], atom[1] + x, atom[2] +y, atom[3] + z]
                    atoms_super.append(new_atom)
    else:
       for i in range(cells[0]):
        for j in range(cells[1]):
            for k in range(cells[2]):
                for atom in  crmg.atoms:
                    new_atom = [atom[0], (atom[1] + i)/float(cells[0]), (atom[2] +j)/float(cells[1]), (atom[3] + k)/float(cells[2])]
                    atoms_super.append(new_atom)

    for i in range(3):
        for j in range(3):
            crmg.cell.latticevectors[i][j] = crmg.cell.latticevectors[i][j] * cells[i] 

    crmg.atoms = atoms_super
