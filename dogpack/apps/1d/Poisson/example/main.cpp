#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "PoissonParams.h"
#include "DogSolver.h"
using namespace std;

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
//  DoGPack
//  -------------------------------
//
//    Lead Developer:  
//             James Rossmanith
//             Iowa State University
//             Department of Mathematics
//             396 Carver Hall
//             Ames, IA 50011
//             rossmani@iastate.edu
// =========================================================================

//---------------------------------------------------------------------------
// Code for solving
//	-u''(x) = f(x) with u'(a) = gamma and u(b) = beta
//---------------------------------------------------------------------------


int main(int argc=0, char**argv=NULL)
{
  // ------------------------------------------------------------
  // Function declarations
  void RunStartScript(string outputdir);
  void ParseArguments(int argc,char**argv,string& outputdir);
  void InitApp(IniDocument& ini_doc);
  void L2Project(int mopt, int istart, int iend,
		 const dTensor2& node,
		 const dTensorBC3& qin, 
		 const dTensorBC3& auxin,  
		 dTensorBC3& Fout,
		 void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
  void ConvertBC(const dTensor2&,const int,const int,const double,
		 const double,const dTensorBC3&,double&,double&);
  void PoissonSolve(const int mstart, const dTensor2& node, const double gamma, 
		    const double beta, const dTensorBC3& Fvals, dTensorBC3& qvals);
  void Output(const dTensor2&,const dTensorBC3&,const dTensorBC3&,double,int,string);
  void SourceTermFunc(const dTensor1& xpts, 
		      const dTensor2& qvals, 
		      const dTensor2& auxvals,
		      dTensor2& fvals);
  // ------------------------------------------------------------
  
  // Output title information
  cout << endl;
  cout << "   ------------------------------------------------   " << endl;
  cout << "   | DoGPack: The Discontinuous Galerkin Package  |   " << endl;
  cout << "   | Developed by the research group of           |   " << endl;
  cout << "   |            James A. Rossmanith               |   " << endl;
  cout << "   |            Department of Mathematics         |   " << endl;
  cout << "   |            University of Wisconsin - Madison |   " << endl;
  cout << "   ------------------------------------------------   " << endl;
  cout << endl;

  // Open parameters.ini file
  ini_doc.initFromFile("parameters.ini");
  IniDocument::Section& ini_sec = ini_doc["dogParams"];
  
  // Get parameters
  dogParams.init();
  dogParamsCart1.init(ini_doc);
  cout << endl;
  
  // Read-in Poisson parameters
  InitApp(ini_doc);
  cout << endl;
  
  // Locally store some parameters
  const int&     nout     = dogParams.get_nout();
  const int*     method   = dogParams.get_method();
  const int      morder   = method[1];
  const int&     meqn     = dogParams.get_meqn();
  const int&     melems   = dogParamsCart1.get_melems();
  const double&  xlow     = dogParamsCart1.get_xlow();
  const double&  xhigh    = dogParamsCart1.get_xhigh();
  const double&  dx       = dogParamsCart1.get_dx();
  const int&     mbc      = dogParamsCart1.get_mbc();
  
  // Parse arguments -- sets directory to which output will be sent,
  //                    the default is "output"
  DogSolver::parse_arguments(argc,argv);
  
  // run startscript
  // to create output directory if it does not exist
  // and copy data files to output directory
  void RunStartScript(int ndims);
  RunStartScript(1);
  
  // Output meqn and nout for plotting purposes
  string qhelp;
  string outputdir = get_outputdir();
  qhelp=outputdir+"/qhelp.dat";
  ofstream out_file(qhelp.c_str(), ios::out);
  
  out_file << setprecision(16);
  out_file << nout << endl << meqn << endl << method[6] << endl;
  out_file << method[1] << endl << melems << endl;
  out_file << setw(24) << scientific << xlow << endl << setw(24)
	   << scientific << xhigh << endl << setw(24) 
	   << scientific << dx << endl;
  out_file.close();
  
  //create needed arrays
  dTensor2    node(melems+1,1);
  dTensorBC3     q(melems,meqn,morder,mbc);
  dTensorBC3 Fvals(melems,1,morder,mbc);
  
  node.set(1,1, xlow );
  for(int n=2; n<=(melems+1); n++) 
    {  node.set(n,1, node.get(n-1,1) + dx );  }
  
  //project f(x) onto the legendre basis
  L2Project(0,1-mbc,melems+mbc,node,Fvals,Fvals,Fvals,&SourceTermFunc);
  
  //Convert input boundary conditions to gamma and beta
  double gamma,beta;
  ConvertBC(node,
	    poissonParams.lft_bc_type,
	    poissonParams.rgt_bc_type,
	    poissonParams.lft_val,
	    poissonParams.rgt_val,
	    Fvals,gamma,beta);
  
  //solve Poisson's equation
  PoissonSolve(1,node,gamma,beta,Fvals,q);
  
  //output poisson solution to file
  Output(node,q,q,0.0,0,outputdir);
  Output(node,q,q,0.0,1,outputdir);

  return 0;
  
}//end of main()
