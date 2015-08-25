#include <fstream>
#include "tensors.h"
using namespace std;

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor1& xpts, 
	       dTensor2& qvals)
{
    int i;
    int numpts=xpts.getsize();
    double x;
    double gamma,mass_ratio;
    double rho_i,u1_i,u2_i,u3_i;
    double P11_i,P12_i,P13_i,P22_i,P23_i,P33_i;
    double rho_e,u1_e,u2_e,u3_e;
    double P11_e,P12_e,P13_e,P22_e,P23_e,P33_e;
    double B1,B2,B3,E1,E2,E3;
    char buffer[256];
    ifstream read_file("param.data", ios::in);

    // Parameters
    read_file >> gamma;
    read_file.getline(buffer,256);
    read_file >> mass_ratio;
    read_file.getline(buffer,256);
    read_file.close();
    
    // Initial conditions
    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);
		
        if(x<0.5e0)
        {
            rho_i = 1.0;
            u1_i  = 0.0;
            u2_i  = 0.0;
            u3_i  = 0.0;
            P11_i = 0.5;
	    P12_i = 0.0;
	    P13_i = 0.0;
	    P22_i = 0.5;
	    P23_i = 0.0;
	    P33_i = 0.5;
	    

	    rho_e = 1.0/mass_ratio;
            u1_e  = 0.0;
            u2_e  = 0.0;
            u3_e  = 0.0;
            P11_e = 0.5;
	    P12_e = 0.0;
	    P13_e = 0.0;
	    P22_e = 0.5;
	    P23_e = 0.0;
	    P33_e = 0.5;

	    B1    = 0.75;
	    B2    = 1.0;
	    B3    = 0.0;
	    E1    = 0.0;
	    E2    = 0.0;
	    E3    = 0.0;	    
        }
        else
        {
            rho_i = 0.125;
            u1_i  = 0.0;
            u2_i  = 0.0;
            u3_i  = 0.0;
	    P11_i = 0.05;
	    P12_i = 0.0;
	    P13_i = 0.0;
	    P22_i = 0.05;
	    P23_i = 0.0;
	    P33_i = 0.05;

	    rho_e = 0.125/mass_ratio;
            u1_e  = 0.0;
            u2_e  = 0.0;
            u3_e  = 0.0;
            P11_e = 0.05;
	    P12_e = 0.0;
	    P13_e = 0.0;
	    P22_e = 0.05;
	    P23_e = 0.0;
	    P33_e = 0.05;

	    B1    = 0.75;
	    B2    =-1.0;
	    B3    = 0.0;
	    E1    = 0.0;
	    E2    = 0.0;
	    E3    = 0.0;
        }

        qvals.set(i,1,  rho_i );       // density (ion)
        qvals.set(i,2,  rho_i*u1_i );  // 1-momentum (ion)
        qvals.set(i,3,  rho_i*u2_i );  // 2-momentum (ion)
        qvals.set(i,4,  rho_i*u3_i );  // 3-momentum (ion)
        qvals.set(i,5,  rho_i*u1_i*u1_i + P11_i );    // energy11 (ion)
	qvals.set(i,6,  rho_i*u1_i*u2_i + P12_i );    // energy12 (ion)
	qvals.set(i,7,  rho_i*u1_i*u3_i + P13_i );    // energy13 (ion)
	qvals.set(i,8,  rho_i*u2_i*u2_i + P22_i );    // energy22 (ion)
	qvals.set(i,9,  rho_i*u2_i*u3_i + P23_i );    // energy23 (ion)
	qvals.set(i,10, rho_i*u3_i*u3_i + P33_i );    // energy33 (ion)
	

	qvals.set(i,11, rho_e );       // density (electron)
        qvals.set(i,12, rho_e*u1_e );  // 1-momentum (electron)
        qvals.set(i,13, rho_e*u2_e );  // 2-momentum (electron)
        qvals.set(i,14, rho_e*u3_e );  // 3-momentum (electron)
	qvals.set(i,15, rho_e*u1_e*u1_e + P11_e );    // energy11 (electron)
	qvals.set(i,16, rho_e*u1_e*u2_e + P12_e );    // energy12 (electron)
	qvals.set(i,17, rho_e*u1_e*u3_e + P13_e );    // energy13 (electron)
	qvals.set(i,18, rho_e*u2_e*u2_e + P22_e );    // energy22 (electron)
	qvals.set(i,19, rho_e*u2_e*u3_e + P23_e );    // energy23 (electron)
	qvals.set(i,20, rho_e*u3_e*u3_e + P33_e );    // energy33 (electron)

	qvals.set(i,21, B1 );         // 1-magnetic field
	qvals.set(i,22, B2 );         // 2-magnetic field
	qvals.set(i,23, B3 );         // 3-magnetic field
	qvals.set(i,24, E1 );         // 1-electric field
	qvals.set(i,25, E2 );         // 2-electric field
	qvals.set(i,26, E3 );         // 3-electric field
		
    }
}
