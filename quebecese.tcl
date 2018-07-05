#just practicing with .tcl scripting and VMD stuff

set x 10
puts "the value of x is:$x"
set text "some text"
puts "the value of text is:$text"

#now we shall try to draw some shapes
#
#NOTE: A molecule or something needs to be loaded already
#in order for us to draw shapes onto it
mol new PRB-0-protein.mae

#NOTE: Objects are given numerical IDs in the order they are drawn
#so the dot would have an ID of 0, the line 1, and the cylinder 3
#(or maybe it start at ID 1? not suuuuuureee.....)

#draw a dot
graphics top point {0 0 10}
#draw a solid line from dot down with a width of 5
graphics top line {-10 0 0} {0 0 0} width 5 style solid
#next object drawn will be color id 3 = orange
graphics top color 3
#drawing cylinders
graphics top cylinder {15 0 0} {15 0 10} radius 10 resolution 69 filled no
#graphics top cylinder {0 0 0} {-15 0 10} radis 5 resolution 60 filled yes


