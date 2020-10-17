load 2.pdb
set defer_builds_mode, 1
load 2.dcd
remove resn HOH
remove resn NA
color silver, all
color orange, resi 381
color cyan, resi 382
show stick, resi 381-382
show sphere, resi 381-382
hide (hydro)
set sphere_transparency, 0.6
set bg_rgb, white
intra_fit all
smooth
set movie_loop, off
