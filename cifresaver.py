""" Read a cif, make a pymatgen structure 
and re-save the resulting structure using symmetrization
The pymatgen project:
https://pymatgen.org/
"""

# avoid CIF export bugs from earlier versions of pymatgen
import pkg_resources
import pymatgen
print("pymatgen version (should be at least 2023.3.23): ",pkg_resources.get_distribution('pymatgen').version)

import os
import sys
import numpy as np

from pymatgen.core import Structure, Lattice, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter

import PySimpleGUI as sg

sg.theme('Dark Blue 3') 

fname = sg.popup_get_file('Input file: Select CIF file to resave',file_types=(("CIF Files", "*.cif"),))
if not fname:
    sg.popup("Cancel", "No filename supplied")
    raise SystemExit("Cancelling: no filename supplied")

filename_in = fname
filename_base, file_extension = os.path.splitext(fname)


filename_out = filename_base +"_resaved"+file_extension

""" output file
"""
layout = [[sg.Text('Enter Output Filename')],
                 [sg.InputText(filename_out), sg.FileBrowse()],
                 [sg.Submit(), sg.Cancel()]]

window = sg.Window('OUTPUT FILE', layout)

event, values = window.read()
window.close()

if values[0]:
    filename_out = values[0] # the first input element is values[0]
else:
    raise SystemExit("Cancelled")




# command line argument is filename
#CIF_FILENAME = sys.argv[1]
CIF_FILENAME = filename_in

# -------- meaning of output filenames -----------------------
filename_noext, ext = os.path.splitext(CIF_FILENAME)

# filename for simply reading and resaving the original CIF
filename_resave = filename_noext+"_resaved" + ext

# structure in P1, plus inverted structure
filename_p1 = filename_noext+"_pmg_P1" + ext
filename_p1_inv = filename_noext+"_P1_INV" + ext

# symmetrized structure from non/inverted coordinates in P1
filename_new = filename_noext+"_symm" + ext
filename_new_inv = filename_noext+"_symm_inv" + ext
# -------------------------------------------------------------


# IMPORT
structure = Structure.from_file(CIF_FILENAME)
print("Structure imported from file: ", CIF_FILENAME)
#print(structure)


# Analyze structure with SpaceGroupAnalyzer 
# BUG??? spg.get_refined_structure() gives P1 structure
# https://pymatgen.org/pymatgen.symmetry.analyzer.html?highlight=spacegroupanalyzer#pymatgen.symmetry.analyzer.SpacegroupAnalyzer
#filename_SGA = filename_noext+"_pmg_SGA" + ext
#spg = SpacegroupAnalyzer(structure)
#spg_info = spg.get_symmetry_dataset()
#print(spg_info)
#periodic_struc = spg.get_refined_structure()
#periodic_struc.to(filename=filename_SGA)


# simply resave the imported structure, hoping that CIF export will produce a valid CIF
#cw = CifWriter(structure, symprec=1.0e-5, significant_figures=6, refine_struct=True)
#cw.write_file(filename_resave)


# now rebuild a new structure (can be modified easier, e.g. inversion, strain)
species = []
coords = []

print("Lattice Parameters:")
lattice_parameters = structure.lattice.parameters
print(lattice_parameters)

for element in structure.composition.elements:
    print(element.symbol, element.Z)

# explicit list of all coordinates for possible manipulation later
print("explicit list of all unit cell atoms:")
atom_count=0
unit_cell_atoms = []
for site in structure.sites:
    occupation = site.species.as_dict()
    site_coords =  site.frac_coords
    for el in site.species:
        coords.append(site_coords.tolist())
        atom_count += 1
        occ = site.species.get_atomic_fraction(el)
        unit_cell_atoms.append([atom_count, el.symbol, el.Z, site_coords, occupation[el.symbol]])
        print(unit_cell_atoms[-1])

#input("Key")

#print("coords:")
#print(coords)

#print("species:")
# https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure
#List of dict of elements/species and occupancies, 
#e.g., [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of disordered structures."
species = []
for idx,site in enumerate(structure.sites):
    occupation = site.species.as_dict()
    for element in occupation:
        species.append({element :  occupation[element] })
    
#print(species)


lattice = Lattice.from_parameters(a=lattice_parameters[0], b=lattice_parameters[1], c=lattice_parameters[2],\
                    alpha=lattice_parameters[3], beta=lattice_parameters[4], gamma=lattice_parameters[5])
                    
structure_new = Structure(lattice, species, coords)
#structure_new.to(filename=filename_p1)
#print(structure_new)

cw = CifWriter(structure_new, symprec=1.0e-5, significant_figures=6, refine_struct=True)
cw.write_file(filename_out)
print("Structure resaved to: ", filename_out)


#structure_inv = Structure(lattice, species, -np.array(coords)+1.0)
#structure_inv.to(filename=filename_p1_inv)
#print(structure_inv)

#cw = CifWriter(structure_inv, symprec=1.0e-5, significant_figures=6, refine_struct=True)
#cw.write_file(filename_new_inv)

input("Press Enter to close...")