import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysisTests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

u = mda.Universe(PSF, DCD)
#u = mda.Universe("/home/bdv1/Desktop/rpn11_ubq.psf", "/home/bdv1/Desktop/rpn11_ubq.dcd")
q1q2 = contacts.q1q2(u, 'name CA', radius=8)
q1q2.run()
plt.ion()

f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
ax[0].legend(loc='best')
ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')



f.show()
plt.pause(20)
