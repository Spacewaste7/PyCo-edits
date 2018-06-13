import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

u = mda.Universe('/home/bdv1/Desktop/rpn11_ubq.psf', '/home/bdv1/Desktop/rpn11_ubq.dcd')

sec_1 = "resname RN11"
section_1 = u.select_atoms(sec_1)
print section_1
sel_1 = "name RN11"
sel_2 = "name UBQ"

selection_1 = u.select_atoms(sel_1)
selection_2 = u.select_atoms(sel_2)

print selection_1
print selection_2

ca1 = contacts.Contacts(u, selection=(sel_1, sel_2), refgroup=(selection_1, selection_2), radius=6.0)

ca1.run()

average_contacts = numpy.mean(ca1.timeseries[:, 1])
print('average contacts = {}'.format(average_contacts))

f, ax = plt.subplots()
ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
ax.set(xlabel='frame', ylabel='fraction of native contacts', title='Native Contacts, average = {:.2f}'.format(average_contacts))
#save('outfile_mdatests')
plt.ion()

f.show()
plt.pause(20)



