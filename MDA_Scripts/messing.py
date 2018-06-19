import MDAnalysis as mda

u = mda.Universe('/home/bdv1/Desktop/rpn11_ubq.psf', '/home/bdv1/Desktop/rpn11_ubq.dcd')


sel1 = u.select_atoms('segid RN11')
print sel1
