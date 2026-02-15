# Hologram_Simulator

This codebase is is for simulating the appearance of a half-parallax hologram, given half-parallax hogel-column images.  It produces frames of an animation of how the resulting hologram would look - accurate within the ray (but not the wave) model of light - by moving the camera around a virtual copy of the hologram and rendering from a few angles.  You can then put the frames together with ffmpeg.

The code is, by my standards, not bad.  I could probably clean it up enough to make it easy to use with a few hours work.  The magic happens inside hologram_simulator.cpp - the main issue is that the various parameters and assumptions are hard coded in there instead of being command line parameters.  For example, right now the line

std::string hogel_folder_name = "Succub_HogelsHR/Hogel_Column_";

means that it is assumed the hogel column images are in a folder called Succub_HogelsHR, and named Hogel_Column_1.png, Hogel_Column_2.png, and so on.  

N_A, N_B hardcode the number of hogels horizontally and vertically, respectively.  Note that it is assumed that there are N_A many hogel column images.  This makes sense - each hogel column image encodes, for a fixed x coordinate, the hogel images with all different possible y coordinates.

DX, DY encode the physical size of the hologram in simulation space.  

N_I is the width of each hogel column image.  The height is necessarily N_B.  I & J are the output size of each frame in the simulation.

You build everything with make -f Makefile_simulator clean, 
make -f Makefile_simulator, and then run with ./hologram_simulator

There are of course implicit assumptions about how the hogel column images were generated.  

