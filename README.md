# sno2_ase_espresso

Very quickly put together example for cassiterite (SnO2).
The directories containg a python script, input JSON file(s) for the script, and a *.png file plotting the result (for some)

Steps:
 - Install Quantum ESPRESSO (QE) on your machine. Linux is straight forward(-ish) but I am new to Mac OS so I added a "dirty" installation guide. (qe_installation_mac.docx)
 - Install ASE with QE-support on your machine. I added a "dirty" installation guide for Mac OS (ase_installation_mac.docx). Linux follows very close....
 - Run the "ecut" example to get a converged basis set via the kinetic energy cutoff
 - Run the "kpts" to get the converged reciprocal-space mesh (for cubic grid of points N x N x N0
 - Run the "surfaces" example to "expose" different faces of the surface
 
 The pseudopotentials are found in the "pseudo" directory
