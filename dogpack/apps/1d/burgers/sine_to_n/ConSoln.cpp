#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart1.h"

using namespace std;

void ConSoln(const int method[], const dTensor2& node, const dTensorBC3& aux,
        const dTensorBC3& q, double t, string outputdir)

{

    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);

    // Grid information:
    const double dx   = dogParamsCart1.get_dx();
    const double xlow = dogParamsCart1.get_xlow();
    const double sqdx = sqrt(dx);

    // Output file names
    string fname1 = outputdir+"/conservation.dat";
    string fname2 = outputdir+"/total-variation.dat";

    ofstream write_file1,write_file2;
    if( t==0 ) 
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
    dTensor2 qsum(meqn, 4); qsum.setall(0.);
    if( dogParams.get_mcapa() < 1 ) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            double maxq = q.get(1,m,1);
            for (int i=1; i<=mx; i++)
            {
                const double x     = xlow + (double(i)-0.5)*dx;
                const double qtmp  = q.get(i,m,1);
                const double abs_q = fabs( q.get(i,m,1) );

                qsum.set(m, 1, qsum.get(m,1) +   dx*qtmp       );  // total mass
                qsum.set(m, 2, qsum.get(m,2) +   dx*abs_q      );  // L1-norm
                qsum.set(m, 3, qsum.get(m,3) + sqdx*qtmp*qtmp  );  // L2-norm
                qsum.set(m, 4, Max( abs_q, qsum.get(m,4) )     );  // L-inf norm

            }
            qsum.set( m, 3, sqrt( qsum.get(m,3) ) );
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            for (int i=1; i<=mx; i++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double qtmp = q.get(i,m,1);
                const double atmp = aux.get(i, 1, dogParams.get_mcapa() );
                const double abs_q = fabs( atmp*q.get(i,m,1) );

                qsum.set(m, 1, qsum.get(m,1) + atmp*dx*qtmp        );  // total mass
                qsum.set(m, 2, qsum.get(m,2) + atmp*dx*abs_q       );  // L1-norm
                qsum.set(m, 3, qsum.get(m,3) + atmp*sqdx*qtmp*qtmp );  // L2-norm
                qsum.set(m, 4, Max( abs_q, qsum.get(m,4) )         );  // L-inf norm
            }
            qsum.set( m, 3, sqrt( qsum.get(m,3) ) );
        }
    }

    // Write data to file
    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file1 << setw(24) << scientific;
        for( int s=1; s <= 4; s++ )
        {
            write_file1 << qsum.get(m,s) << "  ";
        }
    }
    write_file1 << endl;


    // Total variation
    dTensor1 tv(meqn);
    for (int m=1; m<=meqn; m++)
    {
        tv.set(m, fabs(q.get(1,m,1) - q.get(mx,m,1) ) );
        for (int i=2; i<=mx; i++)
        {
            tv.set(m, tv.get(m) + fabs( q.get(i,m,1) - q.get(i-1,m,1) ) );
        }
    }

    // Maximum and minimum value observed during simulation
    dTensor1 MaxQ(meqn); 
    dTensor1 MinQ(meqn);
    for (int m=1; m<=meqn; m++)
    {
        MaxQ.set(m, q.get(1,m,1) );
        MinQ.set(m, q.get(1,m,1) );
        for (int i=2; i<=mx; i++)
        {
            MaxQ.set(m, Max( MaxQ.get(m), q.get(i,m,1) ) );
            MinQ.set(m, Min( MinQ.get(m), q.get(i,m,1) ) );
        }
    }


    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file2 << setw(24) << scientific << tv.get(m) << " ";
        write_file2 << setw(24) << scientific << MinQ.get(m) << " ";
        write_file2 << setw(24) << scientific << MaxQ.get(m) << " ";
    }
    write_file2 << endl;

}
