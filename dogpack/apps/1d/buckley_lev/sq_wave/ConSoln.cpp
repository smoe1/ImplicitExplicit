#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
using namespace std;

void ConSoln(const int method[], const dTensor2& node, const dTensorBC3& aux,
        const dTensorBC3& q, double t, string outputdir)
{
    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1;
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);

    if (t==0) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    if (method[5]<1) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m,0.0);

            for (int i=1; i<=melems; i++)
            {
                double x = node.get(i,1);
                double dtmp = node.get(i+1,1)-node.get(i,1);
                double qtmp = q.get(i,m,1);

                qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
            }
        }
    }
    else // with capacity function
    {
        const int mcapa = method[5];
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

            for (int i=1; i<=melems; i++)
            {
                double x = node.get(i,1);
                double dtmp = node.get(i+1,1)-node.get(i,1);
                double qtmp = q.get(i,m,1);
                double atmp = aux.get(i,mcapa,1);

                qsum.set(m, (qsum.get(m) + atmp*dtmp*qtmp) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;

    string fname2 = outputdir+"/total_variation.dat";
    ofstream write_file2;

    // ----------------------
    // Total Variation
    // ----------------------
    dTensor1 tv(meqn);

    if (t==0) 
    {
        write_file2.open(fname2.c_str(), ofstream::out);
    }
    else
    {
        write_file2.open(fname2.c_str(), ofstream::app);
    }

    const int   kmax = q.getsize(3);
    for (int m=1; m<=meqn; m++)
    {
        tv.set(m,0.0);

        // This loop ignores the boundary cells, and only computes the total
        // variation on the interior of the domain
        for (int i=1; i<=melems-1; i++)
        {

            // evaluate q on the left and right hand side of element:
            double Qr = 0.;
            double Ql = 0.;
            for( int k=1; k <= kmax; k++ )
            {
                Ql += sqrt( 2.0*k-1.0 )*q.get(i,m,k);
                Qr += sqrt( 2.0*k-1.0 )*q.get(i+1,m,k)*pow(-1.0,k-1);
            }
            tv.set(m, tv.get(m) + abs( Ql-Qr ) );
        }
    }

    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file2 << setw(24) << scientific << tv.get(m) << " ";
    }
    write_file2 << endl;

}
