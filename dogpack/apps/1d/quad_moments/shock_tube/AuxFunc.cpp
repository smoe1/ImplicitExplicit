#include <fstream>
#include "tensors.h"
#include "QuadMomentParams.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{
  const int numpts = xpts.getsize();
  const string closure = quadMomentParams.closure;
  
  int cl=0;
  if (closure == "M4")
    {  cl = 1;  }
  else if (closure == "M6")
    {  cl = 2;  }
  else if (closure == "M8")
    {  cl = 3;  }
  
  switch(cl)
    {
    case 1:
      break;
    case 2:
      //////////////////////////// M6 closure //////////////////////////
      for (int i=1;i<=numpts;i++)
        {
	  auxvals.set(i,1,1.0);
	  auxvals.set(i,2,0.5);
	  auxvals.set(i,3,0.25);
	  auxvals.set(i,4,0.5);
	  auxvals.set(i,5,0.25);
	  auxvals.set(i,6,0.25);
        }  
      break;
      
    case 3:
      /////////////////////////// M8 closure ////////////////////////////
      for (int i=1;i<=numpts;i++)
        {
	  auxvals.set(i,1,1.0);
	  auxvals.set(i,2,-1.0);
	  auxvals.set(i,3,0.5);
	  auxvals.set(i,4,0.5);
	  auxvals.set(i,5,0.25);
	  auxvals.set(i,6,0.25);
	  auxvals.set(i,7,0.25);
	  auxvals.set(i,8,0.25);
        }
      break;
      
    }
  
  
}
