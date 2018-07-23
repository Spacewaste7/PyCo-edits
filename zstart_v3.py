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

#acquire distance matrices of backbone atoms for each helix in u for every frame
i = 0
len_frame = len(u.trajectory)
#pos_hel_1 = [] getting matrices of positions was superfolous so I am just going to use these
#loops to acquire the distance matrices
dist_hel_1 = []
for i in range(len_frame):
	u.trajectory[i]
	uatoms = u.select_atoms(hel_1)
	coordU = uatoms.positions
	distU = cdist(coordU, coordU, 'euclidean')
	dist_hel_1.append(distU)
	i = i + 1
i = 0
#pos_hel_2 = []
dist_hel_2 = []
for i in range(len_frame):
        u.trajectory[i]
        uatoms = u.select_atoms(hel_2)
        coordU = uatoms.positions
	distU = cdist(coordU, coordU, 'euclidean')
        dist_hel_2.append(distU)
        i = i + 1
i = 0
#pos_hel_3 = []
dist_hel_3 = []
for i in range(len_frame):
        u.trajectory[i]
        uatoms = u.select_atoms(hel_3)
        coordU = uatoms.positions
	distU = cdist(coordU, coordU, 'euclidean')
        dist_hel_3.append(distU)
        i = i + 1

#acquire distance matrix of backbone atoms for each in helix in ref
distref_hel_1 = []
refatoms1 = ref.select_atoms(hel_1)
coordU = refatoms1.positions
distU = cdist(coordU, coordU, 'euclidean')
distref_hel_1.append(distU)

distref_hel_2 = []
refatoms2 = ref.select_atoms(hel_2)
coordU = refatoms2.positions
distU = cdist(coordU, coordU, 'euclidean')
distref_hel_2.append(distU)

distref_hel_3 = []
refatoms3 = ref.select_atoms(hel_3)
coordU= refatoms3.positions
distU = cdist(coordU, coordU, 'euclidean')
distref_hel_3.append(distU)



#contact determination based on position matrix
#True/False array of contacts for every frame

#make all 0 values (that occur from residues being considered in contact with themselves) into value 50 so they go above the cutoff
dist_hel_1 = numpy.array(dist_hel_1)
dist_hel_2 = numpy.array(dist_hel_2)
dist_hel_3 = numpy.array(dist_hel_3)
distref_hel_1 = numpy.array(distref_hel_1)
distref_hel_2 = numpy.array(distref_hel_2)
distref_hel_3 = numpy.array(distref_hel_3)

dist_hel_1[dist_hel_1 == 0] = 50
dist_hel_2[dist_hel_2 == 0] = 50
dist_hel_3[dist_hel_3 == 0] = 50
distref_hel_1[distref_hel_1 == 0] = 50
distref_hel_2[distref_hel_2 == 0] = 50
distref_hel_3[distref_hel_3 == 0] = 50

#True/False array of contacts/ distance from atoms that are less than 3A
cutoff = 3
tf_hel_1 = dist_hel_1 < cutoff
tf_hel_2 = dist_hel_2 < cutoff
tf_hel_3 = dist_hel_3 < cutoff
tfref_hel_1 = distref_hel_1 < cutoff 
tfref_hel_2 = distref_hel_2 < cutoff
tfref_hel_3 = distref_hel_3 < cutoff

#compare the True/False matrices of trajectory to that of ref structure
print tf_hel_1[0]
print tfref_hel_1[0]
zed = numpy.array(tf_hel_1) == numpy.array(tfref_hel_1)
#print zed[0]


####

#### TRYING TO TAKE POSITIONS OF ALL FALSE VALUES AND THEN MAKE ALL
#### POSITIONS THAT WOULD HAVE HELD FALSE VALUES IN THE ARRAYS INTO FALSES
#### AKA IF THERE WAS A FALSE AT ANY SPOT IN THE COMPARED ARRAYS THEN THIS
#### SCRIPT WILL PUT A FALSE BACK IN THAT POSITION IF IT WAS REMOVED
#(because when you compared them if same positions held a false then
#  it was returned as True)

tfmatrix = zed[0]
tfmatrix = tfmatrix.astype(int)
print('printing base tfmatrix')
tf_hel = tf_hel_1[0].astype(int)
tf_refhel = tfref_hel_1[0].astype(int)

#find x,y coords of 0/false values, then convert list of x,y coords into a
#numpy array to be able to access x,y coords independently of one another
F_coords = zip(*numpy.where(tf_hel==0))
F_coords = numpy.array(F_coords)
REFF_coords = zip(*numpy.where(tf_refhel==0))
REFF_coords = numpy.array(REFF_coords)

#call x and y coords of 0/false values
X = F_coords[:,0]
Xref = REFF_coords[:,0]
Y = F_coords[:,1]
Yref = REFF_coords[:,1]

#change values at positions to 0's/falses
tfmatrix[X,Y] = 0
tfmatrix[Xref,Yref] = 0
print('printing edited tfmatrix')
print tfmatrix

#GREAT, it all works (so faaaaaaar)

####



#uniqueU, countsU = numpy.unique(zed[0], return_counts=True)
#print uniqueU
#print countsU

#comparison to see if they have the same contacts could make True/False
#array for matching Trues

#matrix w/ Trues where both  matrices had T's & F's everywhere else

#read # T's to numeric value

#divide #T's by 2 to account for mirroring in array

#divide #T's by (#T's in ref/2) to get Q

#subtract #T's per frame that match with ref T's value from (#T's in frame/2) 
#to get Z per frame




