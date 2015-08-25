#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "VPDataCart2.h"

void VPDataCart2::init()
{

    // Safe-gaurd in case this is called more than once
    if( is_initialized ){ return; }

    // Dimension various arrays
    const int mx          = dogParamsCart2.get_mx();
    const int my          = dogParamsCart2.get_my();
    const int mbc         = dogParamsCart2.get_mbc();
    const int space_order = dogParams.get_space_order();
    const int time_order  = dogParams.get_time_order();
    const int kmax        = dogParams.get_kmax();
    const int meqn        = dogParams.get_meqn();
    const int kmax1d      = space_order;
    const double       dx = dogParamsCart2.get_dx();
    const double     xlow = dogParamsCart2.get_xlow();
    const double       dy = dogParamsCart2.get_dy();
    const double     ylow = dogParamsCart2.get_ylow();

    // Needed for Poisson solves
    node1d     = new dTensor2(mx+1,1);

    Efield     = new dTensorBC3(mx,1,kmax1d,mbc,1);
    Efield_old = new dTensorBC3(mx,1,kmax1d,mbc,1);
    v1d        = new dTensorBC3(my,1,kmax1d,mbc,1);

    // Kinetic energy
    KE         = new dTensorBC3(mx,1,kmax1d,mbc,1);
    KE_old     = new dTensorBC3(mx,1,kmax1d,mbc,1);

    qold       = new dTensorBC4(mx,my,meqn,kmax,mbc,2);

    for (int i=1; i<=(mx+1); i++)
    {  node1d->set(i,1, xlow + (double(i)-1.0)*dx );  }

//  for (int i=(1-mbc); i<=(mx+mbc); i++)
//  for (int k=1; k<=kmax1d; k++)
//  {
//      Efield->set(i,1,k, 0.0 );
//      Efield_old->set(i,1,k, 0.0 );

//      KE->set(i,1,k, 0.0 );
//      KE_old->set(i,1,k, 0.0 );
//  }
    Efield->setall(0.);
    Efield_old->setall(0.);
    KE->setall(0.);
    KE_old->setall(0.);

    switch( kmax1d )
    {
        case 5:

        case 4:

        case 3:      
            for (int j=(1-mbc); j<=(my+mbc); j++)     
            for (int k=3; k<=kmax1d; k++)
            {
                v1d->set(j,1,k, 0.0 );
            }
        case 2:
            for  (int j=(1-mbc); j<=(my+mbc); j++)
            {
                v1d->set(j,1,2, 0.5*osq3*dy );
            }
        case 1:
            for  (int j=(1-mbc); j<=(my+mbc); j++)
            {
                v1d->set(j,1,1, ylow + (double(j)-0.5)*dy );
            }
            break;

        default:
            printf(" ERROR in VPDataCart2.cpp, kmax1d = %i\n",kmax1d);
            printf("   The current maximum allowed kmax1d is 5\n");
            exit(1);
    }

//  for (int i=(1-mbc); i<=(mx+mbc); i++)
//  for  (int j=(1-mbc); j<=(my+mbc); j++)
//  for (int m=1; m<=meqn; m++)
//  for (int k=1; k<=kmax; k++)
//  {
//      qold->set(i,j,m,k, 0.0 );
//  }
    qold->setall(0.);

    is_initialized = true;

}


VPDataCart2::~VPDataCart2()
{
    if( !is_initialized )
    { return; }

    // free any allocated memory
    delete node1d;
    delete Efield;
    delete Efield_old;
    delete v1d;
    delete KE;
    delete KE_old;
    delete qold;

}
