# solute-diffusion-in-bcc

Program for Solute Diffusion in BCC Solvent

Overview This repository contains a Fortran program developed to analyze the hopping of solutes in a Body-Centered Cubic (BCC) solvent lattice. The program specifically tracks and computes data for four types of hopping events:

Tetra-Tetra (TT) Hopping

Octa-Tetra (OT) Hopping

Tetra-Octa (TO) Hopping

Octa-Octa (OO) Hopping

The program reads in trajectory data, energy values, and void positions to determine the frequency and energy characteristics of these hopping events.

Main Features

Tracks solute movement within a BCC lattice and categorizes the events based on the types of voids between which the solute moves.
Computes and outputs the energy associated with each hopping event, both as calculated by LAMMPS and using the program’s internal calculation methods.
Provides coordination data for each solute during its movement between voids.
Utilizes the Minimum Image Convention for accurate distance calculations in periodic boundary conditions.
Program Structure

The program is structured into the following key modules and subroutines:

Module Parameter: Defines all parameters, variables, and arrays used in the program.
Main Program (BCC): Handles file I/O, data reading, and the main loop for processing trajectory data.
Subroutine voidvoiddistance: Calculates the distance between two void centers.
Subroutine guestvoiddistance: Calculates the distance between a solute and a void center.
Subroutine minimumimageconvention: Adjusts distances using the Minimum Image Convention for periodic boundary conditions.
Subroutine guestenergy: Computes the interaction energy between the solute and surrounding atoms during hopping.
Usage

To run the program, compile the Fortran code with an appropriate compiler (e.g., gfortran) and execute it with the required input files. The program reads trajectory data from input files and outputs the hopping event counts, energy calculations, and coordination data to specified output files.
