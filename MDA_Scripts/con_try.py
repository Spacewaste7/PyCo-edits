# contact analysis
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

u = mda.Universe("/home/bdv1/VMD_SaveFolder/PRB-0-HIE_ProteinOnly.pdb", "/home/bdv1/Desktop/PRB-0-protein-000.dcd")
eitrig = contacts.q1q2(u, 'segid AP1', radius=3)
eitrig.run()

eitrig.save('con_try_contacts.txt')

