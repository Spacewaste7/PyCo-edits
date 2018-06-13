import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy

u = mda.Universe('/home/bdv1/Desktop/rpn11_ubq.psf', '/home/bdv1/Desktop/rpn11_ubq.dcd')

sel_1 = "segid RN11"
sel_2 = "segid UBQ"

selection_1 = u.select_atoms(sel_1)
selection_2 = u.select_atoms(sel_2)
print("Printing selection_1")
print selection_1
print("Printing selection_2")
print selection_2

q1q2 = contacts.q1q2(u, 'segid RN11', radius=8)
q1q2.run()

f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
ax[0].legend(loc='best')
ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')

f.show()
plt.pause(20)
