#imports
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

#specify trajectory and ref universes
u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-025.dcd')
ref = mda.Universe('/home/bdv1/MDA_Scripts/PDB_000-9184.pdb')

#atom selections for helices
hel_1 = "resnum 2:18 and backbone"
hel_2 = "resnum 19:31 and backbone"
hel_3 = "resnum 33:48 and backbone"

#acquire positions matrices of backbone atoms for each helix in u for every frame
i = 0
len_frame = len(u.trajectory)
pos_hel_1 = []
for i in range(len_frame):
	u.trajectory[i]
	uatoms = u.select_atoms(hel_1)
	coordU = uatoms.positions
	pos_hel_1.append(coordU)
	i = i + 1
i = 0
pos_hel_2 = []
for i in range(len_frame):
        u.trajectory[i]
        uatoms = u.select_atoms(hel_2)
        coordU = uatoms.positions
        pos_hel_2.append(coordU)
        i = i + 1
i = 0
pos_hel_3 = []
for i in range(len_frame):
        u.trajectory[i]
        uatoms = u.select_atoms(hel_3)
        coordU = uatoms.positions
        pos_hel_3.append(coordU)
        i = i + 1

#acquire position matrix of backbone atoms for each in helix in ref
refatoms1 = ref.select_atoms(hel_1)
posref_hel_1 = refatoms1.positions

refatoms2 = ref.select_atoms(hel_2)
posref_hel_2 = refatoms2.positions

refatoms3 = ref.select_atoms(hel_3)
posref_hel_3 = refatoms3.positions


#contact determination based on position matrix
#determine distance between every atom and every other atoms per frame


