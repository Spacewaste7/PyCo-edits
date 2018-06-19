import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.tests.datafiles import PSF, DCD
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#u = mda.Universe('/home/bdv1/Desktop/rpn11_ubq.psf', '/home/bdv1/Desktop/rpn11_ubq.dcd')
u = mda.Universe(PSF, DCD)
#sel1 = u.select_atoms('segid RN11')
#sel2 = u.select_atoms('segid UBQ')

sel_basic = "(resname ARG LYS) and (name NH* NZ)"
sel_acidic = "(resname ASP GLU) and (name OE* OD*)"

acidic = u.select_atoms(sel_acidic)
basic = u.select_atoms(sel_basic)

#cmap = distances.distance_array(sel1, sel2, box=sel1.dimensions)
cmap = distances.distance_array(acidic.positions, basic.positions, box=acidic.dimensions)


cutoff = 6
contacts = cmap < cutoff

fig, ax = plt.subplots()
ms = ax.matshow(cmap)
plt.colorbar(ms)
plt.show()


