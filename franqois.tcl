#the frame stuff for VMD
#or something
#

#load the strucuture
mol new PRB-0-protein.mae
#mol load PRB-0-protein-025 $start_dcd 	: not right/working
mol addfile PRB-0-protein-025-EDITS.dcd waitfor all
#the above works properly
set whole [atomselect top all]
set frame 0

#a loop would begin here, apparently
#will see if I need to figure this out later

#add new frame to the trajectory, with coords copied from previous frame
animate dup frame $frame 0
#NOTE: THE ABOVE ONLY WORKS WHEN VMD IS FRESH, E.G. ONLY ONE MOLECULE
#HAS BEEN LOADED AND THUS WOULD HAVE THE MOLID 0, IF IT'S NOT 
#THE FIRST LOADED THEN THE MOLID NUMBER WOULD HAVE TO BE CHANGED FROM
#0 TO WHATEVER IT NEEDS TO BE
#OK
#NEVERMIND
#WHEN COMBINED WITH THE COMMANDS BELOW IT JUST KEEPS ADDING FRAMES TO THE
#MOLECULE WITH MOLID 0	

incr frame
animate goto $frame

#here I need to duplicate a frame from later in PRB-025 to 
#be the very first frame of the trajectory

#end of loop would be here

#write the complete trajectory to disk
#animate write dcd #final_dcd beg 0 end $frame waitfor all

