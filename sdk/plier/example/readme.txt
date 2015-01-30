>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

Copyright (C) 2004 Affymetrix, Inc.

This program is free software; you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation; either version 2 of the License, 
or (at your option) any later version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program;if not, write to the 

Free Software Foundation, Inc., 
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


"Test" is created as a sample program for using the PLIER SDK interfaces and methods.

This file contains a summary of what you will find in each fo the files that
make up the "Test" program.

Test.dsp
     This file (the VC++ project file) contains information at the VC++ project level and
     is used to build the program in VC++.
     
Test.cpp
     This is the main sample source file that defines the PLIER algorithm class (plier_alg) and
     wrapper function (run_plier_wrapper).
     
     The class will read in an input file (input.txt) which contains the input PLIER parameters and
     the input PM and MM data.
     
     The wrapper function (run_plier_wrapper) is used to demonstrate how to wrap up the required PLIER
     interface methods into a single line function call.
     
Makefile
     This is a Linux makefile which is used to compile and link the PLIER SDK and "Test" program.
     


