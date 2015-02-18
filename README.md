# Fortran 95

This program was created by Alexander Cronheim & Aboubakr El Mahdaoui in Fortran. The program simulates an argon gas starting with positions in an fcc lattice and random velocities but with a Maxwell velocity distribution.
Compile with "gfortran argon-box.f90 $(pkg-config --cflags --libs plplotd-f95)". The module plplot is necessary to compile.
