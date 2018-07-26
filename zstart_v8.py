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



#specify trajectory and ref universes
u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-033.dcd')
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
#acquire distance matrices between Helix 1-2, 2-3, and 1-3
dist_hel_12 = []
for i in range(len_frame):
	u.trajectory[i]
	uatoms1 = u.select_atoms(hel_1)
	coordU1 = uatoms1.positions
	uatoms2 = u.select_atoms(hel_2)
	coordU2 = uatoms2.positions	
	distU = cdist(coordU1, coordU2, 'euclidean')
	dist_hel_12.append(distU)
	i = i + 1

dist_hel_23 = []
for i in range(len_frame):
	u.trajectory[i]
	uatoms1 = u.select_atoms(hel_2)
	coordU1 = uatoms1.positions
	uatoms2 = u.select_atoms(hel_3)
	coordU2 = uatoms2.positions
	distU = cdist(coordU1, coordU2, 'euclidean')
	dist_hel_23.append(distU)
	i = i + 1
i = 0
dist_hel_13 = []
for i in range(len_frame):
        u.trajectory[i]
        uatoms1 = u.select_atoms(hel_1)
        coordU1 = uatoms1.positions
	uatoms2 = u.select_atoms(hel_3)
	coordU2 = uatoms2.positions
	distU = cdist(coordU1, coordU2, 'euclidean')
        dist_hel_13.append(distU)
        i = i + 1

#acquire distance matrix of backbone atoms for each in helix in ref in contact with other helixes 
distref_hel_12 = []
refatoms1 = ref.select_atoms(hel_1)
coordU1 = refatoms1.positions
refatoms2 = ref.select_atoms(hel_2)
coordU2 = refatoms2.positions
distU = cdist(coordU1, coordU2, 'euclidean')
distref_hel_12.append(distU)

distref_hel_23 = []
refatoms1 = ref.select_atoms(hel_2)
coordU1 = refatoms1.positions
refatoms2 = ref.select_atoms(hel_3)
coordU2 = refatoms2.positions
distU = cdist(coordU1, coordU2, 'euclidean')
distref_hel_23.append(distU)

distref_hel_13 = []
refatoms1 = ref.select_atoms(hel_1)
coordU1 = refatoms1.positions
refatoms2 = ref.select_atoms(hel_3)
coordU2 = refatoms2.positions
distU = cdist(coordU1, coordU2, 'euclidean')
distref_hel_13.append(distU)

#contact determination based on position matrix
#True/False array of contacts for every frame

#make all 0 values (that occur from residues being considered in contact with themselves) into value 50(randomly chosen large value) so they go above the cutoff
dist_hel_12 = numpy.array(dist_hel_12)
dist_hel_23 = numpy.array(dist_hel_23)
dist_hel_13 = numpy.array(dist_hel_13)
distref_hel_12 = numpy.array(distref_hel_12)
distref_hel_23 = numpy.array(distref_hel_23)
distref_hel_13 = numpy.array(distref_hel_13)

#values being 0 should not really be an issue considering they aren't being compared against themselves but just in case
dist_hel_12[dist_hel_12 == 0] = 50
dist_hel_23[dist_hel_23 == 0] = 50
dist_hel_13[dist_hel_13 == 0] = 50
distref_hel_12[distref_hel_12 == 0] = 50
distref_hel_23[distref_hel_23 == 0] = 50
distref_hel_13[distref_hel_13 == 0] = 50

#True/False array of contacts/ distance from atoms that are less than 6A
cutoff = 6
tf_hel_12 = dist_hel_12 < cutoff
tf_hel_23 = dist_hel_23 < cutoff
tf_hel_13 = dist_hel_13 < cutoff

#uniqueU, countsU = numpy.unique(tf_hel_13[9184], return_counts=True)
#print('printing unqiueU and countsU for 13 9184')
#print uniqueU
#print countsU

tfref_hel_12 = distref_hel_12 < cutoff 
tfref_hel_23 = distref_hel_23 < cutoff
tfref_hel_13 = distref_hel_13 < cutoff

#uniqueU, countsU = numpy.unique(tfref_hel_13, return_counts=True)
#print('printing uniqueU and countsU for ref13')
#print uniqueU 
#print countsU


#compare the True/False matrices of trajectory to that of ref structure
tfmatrix12 = numpy.array(tf_hel_12) == numpy.array(tfref_hel_12)
tfmatrix23 = numpy.array(tf_hel_23) == numpy.array(tfref_hel_23)
tfmatrix13 = numpy.array(tf_hel_13) == numpy.array(tfref_hel_13)

####

#### TRYING TO TAKE POSITIONS OF ALL FALSE VALUES AND THEN MAKE ALL
#### POSITIONS THAT WOULD HAVE HELD FALSE VALUES IN THE ARRAYS INTO FALSES
#### AKA IF THERE WAS A FALSE AT ANY SPOT IN THE COMPARED ARRAYS THEN THIS
#### SCRIPT WILL PUT A FALSE BACK IN THAT POSITION IF IT WAS REMOVED
#(because when you compared them if same positions held a false then
#  it was returned as True)

tfmatrix12 = tfmatrix12.astype(int)
hel_12_compare = tfmatrix12
i = 0
for i in range(len_frame):
	tf_hel = tf_hel_12[i].astype(int)
	tf_refhel = tfref_hel_12.astype(int)
	
	F_coords = zip(*numpy.where(tf_hel==0))
	F_coords = numpy.array(F_coords)
	
	REFF_coords = zip(*numpy.where(tf_refhel==0))
	REFF_coords = numpy.array(REFF_coords)
	
	X = F_coords[:,0]
	Xref = REFF_coords[:,1] #ALRIGHT SO FOR SOME REASON REFF_COORDS HAS AN ADDITONAL FIRST COLUMN FILLED WITH ONLY ZEROS? NOTE: this affects ALL REF. COORDS FOR SOME REASON
	Y = F_coords[:,1]
	Yref = REFF_coords[:,2]

	hel_12_compare[i][X,Y] = 0 #this should work in THEORY
	hel_12_compare[i][Xref,Yref] = 0
	i = i + 1
tfmatrix23 = tfmatrix23.astype(int)
hel_23_compare = tfmatrix23
X = []
Y = []
i = 0
for i in range(len_frame):
        tf_hel = tf_hel_23[i].astype(int)
        tf_refhel = tfref_hel_23.astype(int)

        F_coords = zip(*numpy.where(tf_hel==0))
        F_coords = numpy.array(F_coords)

        REFF_coords = zip(*numpy.where(tf_refhel==0))
        REFF_coords = numpy.array(REFF_coords)
	#print REFF_coords
        X = F_coords[:,0]
        Xref = REFF_coords[:,1]
        Y = F_coords[:,1]
        Yref = REFF_coords[:,2]

        hel_23_compare[i][X,Y] = 0 #this should work in THEORY
        hel_23_compare[i][Xref,Yref] = 0
        i = i + 1
tfmatrix13 = tfmatrix13.astype(int)
hel_13_compare = tfmatrix13
X = []
Y = []
i = 0
for i in range(len_frame):
        tf_hel = tf_hel_13[i].astype(int)
        tf_refhel = tfref_hel_13.astype(int)

        F_coords = zip(*numpy.where(tf_hel==0))
        F_coords = numpy.array(F_coords)

        REFF_coords = zip(*numpy.where(tf_refhel==0))
        REFF_coords = numpy.array(REFF_coords)
	
        X = F_coords[:,0]
        Xref = REFF_coords[:,1]
        Y = F_coords[:,1]
        Yref = REFF_coords[:,2]

        hel_13_compare[i][X,Y] = 0 #this should work in THEORY
        hel_13_compare[i][Xref,Yref] = 0
        i = i + 1
#### 

#read # T's to numeric value
#read for Helix 12 interactions
i = 0
ts_hel12 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(hel_12_compare[i], return_counts=True)
	try:
		ts_hel12.append(countsU[1])
	except:
		ts_hel12.append(0)
	i = i + 1
i = 0
ts_hel23 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(hel_23_compare[i], return_counts=True)
	try:
		ts_hel23.append(countsU[1])
	except:
		ts_hel23.append(0)
	i = i + 1
i = 0
ts_hel13 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(hel_13_compare[i], return_counts=True)
	try:
		ts_hel13.append(countsU[1])
	except:
		ts_hel13.append(0)
	i = i + 1
#print ts_hel13[0]
#print ts_hel13[9184]
#for reference values
uniqueU, countsU = numpy.unique(tfref_hel_12, return_counts=True)
Tsref_hel12 = countsU[1]
uniqueU, countsU = numpy.unique(tfref_hel_23, return_counts=True)
Tsref_hel23 = countsU[1]
uniqueU, countsU = numpy.unique(tfref_hel_13, return_counts=True)
Tsref_hel13 = countsU[1]
#print('printing Tsref_hel13')
#print Tsref_hel13
#t_vals12 = []
#i = 0
#unique12, counts12 = numpy.unique(
#for i in range(len_fram):
#	t_vals12.append(counts12)
#	i = i + 1

#divide #T's by (#T's in ref) to get Q
Q_hel12 = ts_hel12 / Tsref_hel12 
Q_hel23 = ts_hel23 / Tsref_hel23
Q_hel13 = ts_hel13 / Tsref_hel13
 
#subtract #T's per frame that match with ref T's value from (#T's in frame) 
#to get Z per frame
#get 
i = 0
ts_count_hel12 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(tf_hel_12[i], return_counts=True)
	try:
		ts_count_hel12.append(countsU[1])
	except:
		ts_count_hel12.append(0)
	i = i + 1
i = 0
ts_count_hel23 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(tf_hel_23[i], return_counts=True)
	try:
		ts_count_hel23.append(countsU[1])
	except:
		ts_count_hel23.append(0)
	i = i + 1
i = 0
ts_count_hel13 = []
for i in range(len_frame):
	uniqueU, countsU = numpy.unique(tf_hel_13[i], return_counts=True)
	try:
		ts_count_hel13.append(countsU[1])
	except:
		ts_count_hel13.append(0)
	i = i + 1
#get Z per frame
i = 0
Z_hel12 = []
for i in range(len_frame):
	Z = ts_count_hel12[i] - Tsref_hel12
	if Z > 0:
		Z_hel12.append(Z)
	else:
		Z_hel12.append(0)
	i = i + 1
i = 0
Z_hel23 = []
for i in range(len_frame):
	Z = ts_count_hel23[i] - Tsref_hel23
	if Z > 0:	
		Z_hel23.append(Z)
	else:
		Z_hel23.append(0)
	i = i + 1
i = 0
Z_hel13 = []
for i in range(len_frame):
	Z = ts_count_hel13[i] - Tsref_hel13	
	if Z > 0:
		Z_hel13.append(Z)
	else:
		Z_hel13.append(0)
	i = i + 1

#print('printing Q_hel12')
#print Q_hel12

#time for some moving averages oh buddy

N = 50
cumsum, moving_aves12 = [0], []
for i, t in enumerate(Q_hel12, 1):
	cumsum.append(cumsum[i-1] + t)
	if i>=N:
		moving_ave12 = (cumsum[i] - cumsum[i-N])/N
		moving_aves12.append(moving_ave12)

N = 50
cumsum, moving_aves23 = [0], []
for i, t in enumerate(Q_hel23, 1):
	cumsum.append(cumsum[i-1] + t)
	if i>=N:
		moving_ave23 = (cumsum[i] - cumsum[i-N])/N
		moving_aves23.append(moving_ave23)
N = 50
cumsum, moving_aves13 = [0], []
for i, t in enumerate(Q_hel13, 1):
	cumsum.append(cumsum[i-1] + t)
	if i>=N:
		moving_ave13 = (cumsum[i] - cumsum[i-N])/N
		moving_aves13.append(moving_ave13)

#moving average for Zzzzzzz's
N = 50
cumsum, Zmoving_aves12 = [0], []
for i, t in enumerate(Z_hel12, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                Zmoving_ave12 = (cumsum[i] - cumsum[i-N])/N
                Zmoving_aves12.append(Zmoving_ave12)
N = 50
cumsum, Zmoving_aves23 = [0], []
for i, t in enumerate(Z_hel23, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                Zmoving_ave23 = (cumsum[i] - cumsum[i-N])/N
                Zmoving_aves23.append(Zmoving_ave23)
N = 50
cumsum, Zmoving_aves13 = [0], []
for i, t in enumerate(Z_hel13, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                Zmoving_ave13 = (cumsum[i] - cumsum[i-N])/N
                Zmoving_aves13.append(Zmoving_ave13)






xaxis = [i for i in range(10000)]
f, ax = plt.subplots()

#graph of helix interactions with moving average
xaxis50 = [i for i in range(9951)]
#ax.plot(xaxis50, moving_aves12, 'b', label1='Helix 1-2 interactions', xaxis50, moving_aves23, 'r', xaxis50, moving_aves13, 'g')
ax.plot(xaxis50, moving_aves12, 'b', label='Helix 1-2 interactions')
ax.plot(xaxis50, moving_aves23, 'g', label='Helix 2-3 interactions')
ax.plot(xaxis50, moving_aves13, 'r', label='Helix 1-3 interactions')

#now going to try making a histogram of Z-values
#####	via this testing it has come to my attention I don't know how histograms work
#x = num_bins = 100
#n, bins, patches = plt.hist(Z_hel12, num_bins, facecolor='blue', alpha=0.5)
#plt.show()
#plt.pause(240)

#ax.plot(xaxis50, Zmoving_aves12, 'b', xaxis50, Zmoving_aves23, 'r', xaxis50, Zmoving_aves13, 'g')
#ax.plot(xaxis, Z_hel12, 'b', xaxis, Z_hel23, 'r', xaxis, Z_hel13, 'g')
#hel12_patch = mpatches.Patch(color='blue', label='Helix 1-2 interactions')
#plt.legend(handles=[hel12_patch])
#hel23_patch = mpatches.Patch(color='green', label='Helix 2-3 interactions')
#plt.legend(handles=[hel23_patch])
#hel13_patch = mpatches.Patch(color='red', label='Helix 1-3 interactions')
#plt.legend(handles=[hel13_patch])
#plt.legend([hel12_patch, 
#ax.plt.legend((line1, line2, line3), ('label1','label2','label3')

ax.legend()
ax.set(xlabel='frame', ylabel='Q, fraction of native contacts', title='PRB-033 Fraction of Native Contacts between helices')
ax.set_ylim(0, 1.0)
f.show()
plt.pause(240)










