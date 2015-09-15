#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
        dTensor2& auxvals)
{

    double xmin=100.0;
    double xmax=-100.0;

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        xmin=Min(xmin,x);
        xmax=Max(xmax,x);
        {auxvals.set(i,1, 1.0e0 );}
        {auxvals.set(i,2, 0.0 );}
        {auxvals.set(i,3, -1.0 );}
        //if(x>=0.5 && x<0.7+sqrt(2.0)/10.0)
        //if(x>=0.5 && x<0.7)
        //{auxvals.set(i,1,1.0);}
    }

    for (int i=1; i<=numpts; i++)
    {
        {auxvals.set(i,1, 1.0e0 );}
        {auxvals.set(i,2, 0.0 );}
        //if(x>=0.5 && x<0.7+sqrt(2.0)/10.0)
	if(xmin>=0.5)//&& xmax<0.7)
        {auxvals.set(i,1,1.0);}
        //{auxvals.set(i,1,-0.5);}
    }




    //printf("xmax=%e and xmin=%e %d \n",xmin,xmax,(xmin<0.5 && xmax>0.5) || (xmin<0.7 && xmax>0.7));

    
    if((xmin<0.5 && xmax>0.5))
    {
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        {auxvals.set(i,2, 1.0 );}
        {auxvals.set(i,3, 0.5 );}
    }
    }  

    /*
    if((xmin<0.7 && xmax>0.7))
    {
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        {auxvals.set(i,2, 1.0 );}
        {auxvals.set(i,3, 0.7 );}
    }
    }*/

    
}
