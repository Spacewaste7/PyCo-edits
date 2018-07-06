#load in topology of PRB
mol new PRB-0-protein.mae
#load in trajectory of 000
mol addfile PRB-0-protein-000.dcd waitfor all

#load in topology of PRB
#mol new PRB-0-protein.mae	:	just working with 000 for now
#load in trajectory of 025
#mol addfile PRB-0-protein-025-EDITS.dcd waitfor all	:	^

#setting variables
set whole [atomselect top all]
#so we want 'frame' to be the folded structure
#frame 9184 of trajectory 000 looks good
set frame 9184

#duplicate frame 9184 of molID 0(for 000)
#animate dup frame $frame 0 
#alright so the above just moved the frame to 
#the beggining of 000
#
#let's just work with 000 for now

#testing adding more frames
#animate dup frame $frame 0
#animate dup frame $frame 0
#animate dup frame $frame 0
#succesfuly added frames to end of 000

#will now see if changing the frame vmd is set at will 
#change where the frames are duplicated too
animate goto 50
animate dup frame $frame 0
animate dup frame $frame 0
animate dup frame $frame 0
#alright so that doesn't work
animate dup 50 $frame 0


animate write dcd EDITED-PRB-0.dcd beg 0 waitfor all
