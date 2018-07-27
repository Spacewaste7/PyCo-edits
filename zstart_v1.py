#This script is for getting proper fraction of native contact values as
#well as Z values, the value of non-native contacts


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


#distance matrix comprised of backbone atoms for residues
#SO here is where we use MDA to create distance matrix for all of our selections
#will read from Universes 
U_hel_1 = u.select_atoms(hel_1)
coordU_hel_1 = U_hel_1.positions
coordU = u.coord.positions
coordREF = ref.coord.positions
print coordU_hel_1
print coordU


#contact determination based on distance matrix
#determine distance between every atom and every other atom PER FRAME
#going to need a loop
distU = cdist(coordU, coordU, 'euclidean') #have you guys heard about LOOPER?
#print distU
#distance matrix for reference structure
distREF = cdist(coordREF, coordREF, 'euclidean')
#print distREF
#contact analysis
#so make all Y array values below a cutoff into a true
cutoff = 3
cont_arrU = distU < cutoff
cont_arrREF = distREF < cutoff
#print cont_arrU
#print cont_arrREF


#with True/False array, get # Trues in trajectory per frame
uniqueU, countsU = numpy.unique(cont_arrU, return_counts=True)
print uniqueU[1]
print countsU[1]
#print dict(zip(uniqueU, countsU))
#get #T's in ref. struct
#uniqueREF, countsREF = numpy.unique()
#print dict(zip(uniqueREF, countsREF))


#compare each frames T/F array to that of the ref. structure & count
#how many T's are the same, divide that num by # T's in ref. struct.
#to get Q for that frame


#difference in # T's per frame & # matching T's per frame = Z value


#graphing of multiple things


