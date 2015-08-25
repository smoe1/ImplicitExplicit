#include "tensors.h"
#include <cmath>
#include "QuadMomentParams.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, 
        dTensor2& qvals)
{
    const int numpts = xpts.getsize();
    const string closure = quadMomentParams.closure;

    int cl = 0;
    if (closure == "M4")
    {  cl = 1;  }
    else if (closure == "M6")
    {  cl = 2;  }
    else if (closure == "M8")
    {  cl = 3;  }

    // Initial conditions
    switch(cl)
    {
        case 1:

            ////////////////////////// M4 closure/////////////////////////////  
            for (int i=1; i<=numpts; i++)
            {
                double rad = xpts.get(i);
                double rho,press,u,energy,q,m3;

                if(rad<0.5e0)
                {
                    rho   =  1.0;
                    u     =  0.0;
                    press =  1.0;
                    q     =  0.0;
                }
                else
                {
                    rho   =  0.125;//3.0;
                    u     =  0.0;
                    press =  0.1;//3.0;
                    q     =  0.0;
                }

                energy = press + rho*pow(u,2);
                m3 = q + rho*pow(u,3) + 3.0*press*u;

                qvals.set(i,1, rho );       // density
                qvals.set(i,2, rho*u );     // momentum
                qvals.set(i,3, energy );    // energy
                qvals.set(i,4, m3);         // 3rd-order moment
            }
            break;

        case 2:

            ////////////////////////// M6 closure/////////////////////////////  
            for (int i=1; i<=numpts; i++)
            {
                double rad = xpts.get(i);
                double rho,press,u,energy,q,m3,r,m4,s,m5;

                if(rad>0.5e0)
                {
                    rho    =  1.0;
                    u      =  0.0;
                    press  =  1.0;
                    q      =  0.0;
                    r      =  3.0;
                    s      =  0.0;
                }
                else 
                { 
                    rho    =  2.0;
                    u      =  0.0;
                    press  =  2.0;
                    q      =  0.0;
                    r      =  6.0;
                    s      =  0.0;
                }

                energy = press + rho*pow(u,2);
                m3 = q + rho*pow(u,3) + 3.0*press*u;
                m4 = r + 4.0*q*u + 6.0*press*pow(u,2) + rho*pow(u,4);
                m5 = s + 5.0*r*u + 10.0*q*pow(u,2) + 10.0*press*pow(u,3) + rho*pow(u,5);

                qvals.set(i,1, rho );       // density
                qvals.set(i,2, rho*u );     // momentum
                qvals.set(i,3, energy );    // energy
                qvals.set(i,4, m3);         // 3rd-order moment
                qvals.set(i,5, m4);         // 4th-order moment
                qvals.set(i,6, m5);         // 5th-order 
            }
            break;

        case 3:

            //////////////////////////////// M8 closure /////////////////////////////////
            for (int i=1; i<=numpts; i++)
            {
                double rad = xpts.get(i);
                double rho,press,u,energy,q,m3,r,m4,s,m5,delta,m6,tau,m7;

                if(rad>0.5e0)
                {
                    rho    =  1.0;
                    u      =  0.0;
                    press  =  1.0;
                    q      =  0.0;
                    r      =  2.0;
                    s      =  0.0;
                    delta  =  12.0;
                    tau    =  0.0;
                }
                else 
                { 
                    rho    =  1.2;
                    u      =  0.0;
                    press  =  1.2;
                    q      =  0.0;
                    r      =  2.4;
                    s      =  0.0;
                    delta  =  14.4;
                    tau    =  0.0;
                }

                energy = press+rho*pow(u,2);
                m3 = q+rho*pow(u,3)+3.0*press*u;
                m4 = r+4.0*q*u+6.0*press*pow(u,2)+rho*pow(u,4);
                m5 = s+5.0*r*u+10.0*q*pow(u,2)+10.0*press*pow(u,3)+rho*pow(u,5);
                m6 = delta+6.0*s*u+15.0*r*pow(u,2)+20.0*q*pow(u,3)+15*press*pow(u,4)+rho*pow(u,6);
                m7 = tau+7.0*delta*u+21.0*s*pow(u,2)+35.0*r*pow(u,3)+35.0*q*pow(u,4)+21.0*press*pow(u,5)+rho*pow(u,7);


                qvals.set(i,1, rho );       // density
                qvals.set(i,2, rho*u );     // momentum
                qvals.set(i,3, energy );    // energy
                qvals.set(i,4, m3);         // 3rd-order moment
                qvals.set(i,5, m4);         // 4th-order moment
                qvals.set(i,6, m5);         // 5th-order 
                qvals.set(i,7, m6);         // 6th-order
                qvals.set(i,8, m7);         // 7th-order
            }
            break;

    }

}
