import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD

import scipy
from scipy.interpolate import spline

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-025.dcd')

ref = mda.Universe('/home/bdv1/MDA_Scripts/PDB_000-9184.pdb')

hel_1 = "resnum 2:18"
hel_2 = "resnum 20:31"
hel_3 = "resnum 33:47"

helix_1 = ref.select_atoms(hel_1)
helix_2 = ref.select_atoms(hel_2)
helix_3 = ref.select_atoms(hel_3)

ca12 = contacts.Contacts(u, selection=(hel_1, hel_2), refgroup=(helix_1, helix_2), radius=2000, method=contacts.radius_cut_q, kwargs={"radius":15.0}, start=4200, stop=8200).run()
ca23 = contacts.Contacts(u, selection=(hel_2, hel_3), refgroup=(helix_2, helix_3), radius=2000, method=contacts.radius_cut_q, kwargs={"radius":15.0}, start=4200, stop=8200).run()
ca13 = contacts.Contacts(u, selection=(hel_1, hel_3), refgroup=(helix_1, helix_3), radius=2000, method=contacts.radius_cut_q, kwargs={"radius":15.0}, start=4200, stop=8200).run()


x12 = ca12.timeseries[:,0]
y12 = ca12.timeseries[:,1]

x23 = ca23.timeseries[:,0]
y23 = ca23.timeseries[:,1]

x13 = ca13.timeseries[:,0]
y13 = ca13.timeseries[:,1]

N = 50
cumsum, moving_aves12 = [0], []
for i, t in enumerate(y12, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                moving_ave12 = (cumsum[i] - cumsum[i-N])/N
                moving_aves12.append(moving_ave12)

N = 50
cumsum, moving_aves23 = [0], []
for i, t in enumerate(y23, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                moving_ave23 = (cumsum[i] - cumsum[i-N])/N
                moving_aves23.append(moving_ave23)

N = 50
cumsum, moving_aves13 = [0], []
for i, t in enumerate(y13, 1):
        cumsum.append(cumsum[i-1] + t)
        if i>=N:
                moving_ave13 = (cumsum[i] - cumsum[i-N])/N
                moving_aves13.append(moving_ave13)


#print cumsum
#print moving_aves

f, ax = plt.subplots()
#ax.plot(ca12.timeseries[1:320, 0], ca12.timeseries[1:320, 1], 'b')

#ax.plot(x_smooth, y_smooth, 'r')
ax.plot(x13[1:3952], moving_aves12, 'b', x23[1:3952], moving_aves23, 'r', x13[1:3952], moving_aves13, 'g')
ax.set(xlabel='frame', ylabel='fraction of native contacts', title='PRB-025 Native Contacts, frames 4200-8200')
ax.set_ylim(0, 1.0)
f.show()
plt.pause(30)

