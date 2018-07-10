import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-025.dcd')

#ref = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-000.dcd')
#ref.trajectory[9184]
#Ah, make sure that ref uses frame 9184 of traj 000,not another one, that would be incorrect

ref = mda.Universe('/home/bdv1/MDA_Scripts/PDB_000-9184.pdb')
#reference structure created using MDA write command
#first, testing to see if the PDB file of default results in same
#ref = mda.Universe('/home/bdv1/MDA_Scripts/PDB_000-default.pdb')
#it does, atleast at radius 20

#ref = mda.Universe('/home/bdv1/MDA_Scripts/VMD_test_newdrop/PRB-0-protein_autopsf_formatted(good one?).pdb')
# wait why would this work I just made this I have to make it into the one
# frame thing again what the heck
hel_1 = "resnum 2:18"
hel_2 = "resnum 20:31"
hel_3 = "resnum 33:47"

helix_1 = ref.select_atoms(hel_1)
helix_2 = ref.select_atoms(hel_2)
helix_3 = ref.select_atoms(hel_3)

ca1 = contacts.Contacts(u, selection=(hel_1, hel_2), refgroup=(helix_1, helix_2), radius=2000, method=contacts.radius_cut_q, kwargs={"radius":15.0}, stop=9999).run()

print ca1.timeseries[:, 0]
print ca1.timeseries[:, 1]

f, ax = plt.subplots()
ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
#ax.set(xlabel='frame', ylable='fraction of native contacts', title='PRB-025 Native Contacts')
ax.set_ylim(0, 1.0)
f.show()
plt.pause(30)




