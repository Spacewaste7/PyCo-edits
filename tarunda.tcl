mol new PRB-0-protein.mae
mol addfile PRB-0-protein-000.dcd waitfor all

set sel0 [atomselect 0 all]
puts "hon"
puts $sel0

#parray sel0

