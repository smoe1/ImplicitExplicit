#include "PlasmaSolverCart2.h"
#include "IniDocument.h"
// =========================================================================
//
//  Copyright J.A. Rossmanith
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
//  -------------------------------
//  DoGPack Version 1.0,  Jan. 2007
//  -------------------------------
//
//    Author:  James A. Rossmanith
//             Department of Mathematics
//             University of Wisconsin
//             Madison, WI 53706
//             rossmani@math.wisc.edu
// =========================================================================

int main_plasma(int argc, char* argv[])
{
  //
  // NOTE: You should not have to modify this part of the code.
  //       To change parameters, modify the following files:
  //            1. parameters.ini -- basic data file, can modify
  //                  # of grid points, time step, order of
  //                  accuracy in both space and time, etc...
  //            2. QinitFunc.cpp -- initial condition file
  //            3. AuxFunc.cpp -- auxiliary variable file
  //            4. SourceTermFunc.pp -- source term file
  //            5. FluxFunc.cpp -- flux function file
  //            6. SetWaveSpd.cpp -- eigenvalues of flux Jacobian file
  //            7. ProjectLeftEig.cpp -- left eigenvectors of flux Jacobian
  //            8. ProjectLeftEig.cpp -- right eigenvectors of flux Jacobian
  //            9. SetBndValues.cpp -- boundary conditions files
  //
  
  ini_doc.initFromFile("parameters.ini");
  IniDocument::Section& ini_sec = ini_doc["dogParams"];

  DogSolver::parse_arguments(argc,argv);
  void RunStartScript(int);
  RunStartScript(2);

  PlasmaSolverCart2 plasmaSolver;
  int m = plasmaSolver.RunDogpack();

  return m;
}
