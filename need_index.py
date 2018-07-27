#imports
from __future__ import division
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
#from MDAnalysis.tests.datafiles import PSF, DCD
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy
#from __future__ import division
import scipy
from scipy.interpolate import spline

#trying blockwise_view
from blockwise_view import blockwise_view

u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-033.dcd')
ref = mda.Universe('/home/bdv1/MDA_Scripts/PDB_000-9184.pdb')

hel_1 = "resnum 2:18 and backbone"
hel_2 = "resnum 19:31 and backbone"
hel_3 = "resnum 33:48 and backbone"


i = 0
len_frame = len(u.trajectory)

distref_hel_12 = []
refatoms1 = ref.select_atoms(hel_1)
coordU1 = refatoms1.positions
refatoms2 = ref.select_atoms(hel_2)
coordU2 = refatoms2.positions
distU = cdist(coordU1, coordU2, 'euclidean')
distref_hel_12.append(distU)

#print('printing coordU1')
#print coordU1
#print refatoms1.resnums
#print('printing coordU2')
#print coordU2
#print refatoms2.resnums
#print('printing coordU3')
refatoms3 = ref.select_atoms(hel_3)
coordU3 = refatoms3.positions
#print coordU3
#print refatoms3.resnums

distref_hel_12 = numpy.array(distref_hel_12)

distref_hel_12[distref_hel_12 == 0] = 50

cutoff = 6

tfref_hel_12 = distref_hel_12 < cutoff

#print(tfref_hel_12)

#creating a def to cut block into smaller bits

#trying things out
#c = tfref_hel_12
#print(c)
#blockshaped(TF_array, nresidues_2ndhel, nresidues_1sthel)
#print(blockshaped(c, 12, 16))

#a = numpy.arange(1,21).reshape(4,5)
#print a
#blocks = blockwise_view(a, blockshape=(2,2), require_aligned_blocks=False)
#print blocks


#print('DOIN THE BLOCKWISE THING HERE WE GO')
#c = tfref_hel_12
#print c
#tfblocks = blockwise_view(c, blockshape=(1,1), require_aligned_blocks=True)
#print tfblocks
#well that didn't work, maybe try givin them a non True False array
#niet

a = numpy.arange(1,17).reshape(4,4)
print a
#x = numpy.arange(8.0)
#print x
n = 4
slice_list = [slice(k, l) for k in range(0, n) for l in range(k, n)]

results = [a[sl, sl] for sl in slice_list]
print results


	
