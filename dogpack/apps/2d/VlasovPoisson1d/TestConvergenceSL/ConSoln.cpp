#include<iomanip>
#include<fstream>
#include<iostream>
#include<cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
const char* get_outputdir();

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Function for computing conserved quantities for the Vlasov-Poisson equation.
//
//      e(:,1) = abs(f)
//      e(:,2) = f^2
//      e(:,3) = v^2 * f // STILL NEED TO ADD IN E FIELD to get energy //
//      e(:,4) = abs(f) * log( f )
///////////////////////////////////////////////////////////////////////////////
void ConservedFunc(const dTensor2& xpts, const dTensor2& qvals, 
        const dTensor2& auxvals, dTensor2& e)
{

    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double v = xpts.get(i,2);

        double f = qvals.get(i,1);

        e.set(i, 1, fabs( f )       );
        e.set(i, 2, pow(f,2)        );
        e.set(i, 3, 0.5*v*v*f       );

//assert( f > -1e-10 );

        // safety check for taking the log here ... 
        if( fabs(f) < 1e-15 )
        { e.set(i, 4, 0.0 ); }
        else
        { e.set(i, 4, -fabs( f ) * log( fabs(f) ) ); }

    }

}

void ConSoln(const dTensorBC4& aux, const dTensorBC4& q, double t)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    string fname1 = string(get_outputdir())+"/conservation.dat";
    string fname2 = string(get_outputdir())+"/Efield.dat";

    ofstream write_file1,write_file2;

    dTensor1 qsum(meqn);

    if (t==0) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
        write_file2.open(fname2.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
        write_file2.open(fname2.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    if (dogParams.get_mcapa()<1) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m,0.0);

            for (int i=1; i<=mx; i++)	    
            for (int j=1; j<=my; j++)
            {
                qsum.set(m, qsum.get(m) + dx*dy*q.get(i,j,m,1) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (fabs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;

    write_file1.close();

    ///////////////////////////////////////////////////////////////////////////
    // Conserved Vlasov-Poisson Quantities                                   //
    //        ||f||_1, ||f||_2, Energy, Entropy                              //
    ///////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // Electric Field 
    const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d  = int(sqrt(mpoints));
    const int mcons = 4;

    dTensorBC3 Evals(mx, meqn, kmax1d, mbc,1);
    void ComputeElecField(double t, const dTensor2& node1d,
            const dTensorBC4& qvals, dTensorBC3& Evals);

    // save 1d grid points
    dTensor2 node1d(mx+1, 1);
    for(int i=1; i<=(mx+1); i++)
    { node1d.set(i, 1, dogParamsCart2.get_xl(i)); }

    ComputeElecField(t, node1d, q, Evals);
    //////////////////////////////////////////////////////////////////////////

    ///////////// TODO -- THIS IS A WASTE HERE, BECAUSE WE ONLY NEED THE FIRST
    /////////////         WEIGHT!
    void L2Project(const int istart, 
		   const int iend, 
		   const int jstart, 
		   const int jend,
		   const int QuadOrder, 
		   const int BasisOrder_qin,
		   const int BasisOrder_auxin,
		   const int BasisOrder_fout,		  
		   const dTensorBC4* qin,
		   const dTensorBC4* auxin, 
		   dTensorBC4* fout,
		   void (*Func)(const dTensor2&,const dTensor2&,
				const dTensor2&,dTensor2&));
    dTensorBC4 ConsData2d(mx, my, mcons, kmax,mbc); 
    const int space_order = dogParams.get_space_order();
    L2Project(1, mx, 1, my, space_order, space_order, space_order,
	      space_order, &q, &aux, &ConsData2d, &ConservedFunc );

    // electric field, entropy and energy
    dTensor1 ConsData(mcons); 
    ConsData.setall(0.);

    // Compute the integral of all the 2d quantities:
    for (int m=1; m<=mcons; m++)
    {
        for (int i=1; i<=mx; i++)	    
        for (int j=1; j<=my; j++)
        {
            ConsData.set(m, ConsData.get(m) + dx*dy*ConsData2d.get(i,j,m,1) );
        }
    }

    // Compute the integral of all the 1d quanties:
    //        \int E^2 \ dx
    double E2 = 0.0;
    for(int i=1; i <= mx; i++ )
    for(int k=1; k <= kmax1d; k++ )
    { E2 += dx * pow( Evals.get(i,1,k), 2 ); }

    // Add in E to the total Energy:
    ConsData.set(3, ConsData.get(3) + 0.5 * E2 );

    // print time:
    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";

    // print ||E||_2 first:
    if (fabs(E2) < 1.0e-99) { E2 = 0.0; }
    write_file2 << setw(24) << scientific << sqrt( E2 ) << " ";

    // print all the other conserved quantities:
    for (int m=1; m<=mcons; m++)
    {
        if (fabs(ConsData.get(m)) < 1.0e-99) {ConsData.set(m, 0.0);}
        write_file2 << setw(24) << scientific << ConsData.get(m) << " ";
    }
    write_file2 << endl;

    write_file2.close();

}
