import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

u = mda.Universe('/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb', '/home/bdv1/Desktop/PRB-0-protein-000.dcd')

sel_1 = "segid AP1"
sel_2 = "segid AP2"

selection_1 = u.select_atoms(sel_1)
selection_2 = u.select_atoms(sel_2)
#print selection_1
#print selection_2

#q1q2 = contacts.q1q2(u, selection=(sel_1, sel_2), radius=6)
#q1q2.run()

#ca1 = contacts.Contacts(u, selection=(sel_1, sel_2), refgroup=(selection_1, selection_2), radius=6.0)
ca1 = contacts.Contacts(u, selection=(sel_1, sel_2), refgroup=(selection_1, selection_2), radius=6.0, start=0, stop=1001)

ca1.run()

#f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
#ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
#ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
#ax[0].legend(loc='best')
#ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')

average_contacts = numpy.mean(ca1.timeseries[:, 1])
print('average contacts = {}'.format(average_contacts))


f, ax = plt.subplots()
ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
ax.set(xlabel='frame', ylabel='fraction of native contacts',
       title='Native Contacts, average = {:.2f}'.format(average_contacts))
# HAVE matplotlib.use('PDF') then have matplotlib save the figure after it's been made to a pdf then I won;t have to worry about the creation of 
f.show()
#plt.pause(20)
