#load first structure
mol new PRB-0-protein.mae
mol addfile PRB-0-protein-000.dcd waitfor all
#load second structure
#the one which will be edited
#mol new PRB-0-protein.mae
#mol addfile PRB-0-protein-025-EDITS.dcd waitfor all

#load a SECOND 000 structure and trajectory
#the one which will be edited
mol new PRB-0-protein.mae
mol addfile PRB-0-protein-000.dcd waitfor all

set sel0 [atomselect 0 all]
set sel1 [atomselect 1 all]

#select ref structure frame to input to trajectory
$sel0 frame 9184
#select frame to place reference frame
$sel1 frame 0
#set coords of frame to those of ref structure frame
$sel1 set {x y z} [$sel0 get {x y z}]

#delete first frame with unfolded ref structure
#animate delete beg 0 end  0 skip -1 1	: NOT NEEDED

#make new .dcd file
animate write dcd EDITED-PRB-0.dcd beg 0 waitfor all 1 

#mol new PRB-0-protein.mae
#mol addfile EDITED-PRB-0.dcd waitfor all
