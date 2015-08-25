///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

void FiniteDiff(int mestart, int meend, const dTensorBC3& q, dTensorBC3& qout);
void SetPeriodicBndValues1D( dTensorBC3& q);
void ComputeElecField(double t, const dTensor2& node1d,
    const dTensorBC4& qvals, dTensorBC3& Evals);
void ConvertQ2dToQ1d(const int &mopt, int istart, int iend, 
                     const dTensorBC4& qin, dTensorBC3& qout);
void IntegrateQ1dMoment1(const dTensorBC4& q2d, dTensorBC3& q1d);
void IntegrateQ1dMoment2(const dTensorBC4& q2d, dTensorBC3& q1d);
void IntegrateQ1d(const int mopt, const dTensorBC4& q2d, dTensorBC3& q1d);
void MultiplyFunctions(const dTensorBC3& Q1, const dTensorBC3& Q2,
        dTensorBC3& qnew);
void IntegrateQ1dMoment3(const dTensorBC4& q2d, dTensorBC3& q1d);

// The purpose of this routine is to compute the necessary terms for the
// expansion in 
//
//    \bar{E}(t) = E^n + (t-tn) E_t + (t-tn)^2 / 2 * E_tt + (t-tn)^3 / 3! E_ttt
//
// This involves computing moments and derivatives of the distribution function,
// and is called once per time step.  This routine need only be called when we
// wish to perform a method that is 4th order accurate in time, because Strang
// splitting gets 2nd order accuracy for free.
//
// Results are saved into sl_state->aux1d.  See SLState.h for what that
// auxiliary array represents
void InitHybridState( 
    const dTensorBC4& q, const dTensorBC4& aux, SL_state& sl_state )
{

    // The hybrid method requires running this routine EVERY time we initialize
    // the velocities ...

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    const int kmax1d  = dogParams.get_space_order();
    const int mpoints = kmax1d*kmax1d;

    const double tn     = sl_state.tn;

    //////////////////////// Compute Electric Field E(t) ////////////////////
    // save 1d grid points ( this is used for poisson solve )
    dTensorBC3 Enow(mx, meqn, kmax1d, mbc, 1);
    ComputeElecField( tn, *sl_state.node1d, *sl_state.qnew, Enow);

    //////////// Necessary terms for Et = -M1 + g(x,t) /////////////////////
    dTensorBC3 M1(mx, meqn, kmax1d,mbc,1); 
    IntegrateQ1dMoment1(q, M1);  // first moment

    //////////// Necessary terms for Ett = 2*KE_x - rho*E + g_t - psi_u ///////
    dTensorBC3 KE(mx, meqn, kmax1d,mbc,1); 
    IntegrateQ1dMoment2(q, KE);  // 1/2 * \int v^2 * f dv //

    SetPeriodicBndValues1D( KE );
    SetPeriodicBndValues1D( Enow );
    SetPeriodicBndValues1D( M1 );

    // do something to compute KE_x ... //
    dTensorBC3 gradKE(mx,2,kmax1d,mbc,1);
    FiniteDiff( 1, 2, KE, gradKE ); // compute KE_x and KE_xx //

    dTensorBC3 rho(mx,meqn,kmax1d,mbc,1);
    dTensorBC3 prod(mx,meqn,kmax1d,mbc,1);
    IntegrateQ1d( 1, q, rho );
    MultiplyFunctions( rho, Enow, prod );

    /////////////////////// Save E, Et, Ett, //////////////////////////////////
    for( int i = 1; i <= mx; i++ )
    for( int k = 1; k <= kmax1d; k ++ )
    {
        double E     =  Enow.get  ( i, 1, k );
        double Et    = -M1.get ( i, 1, k );  // + g(x,t)
        double Ett   = -prod.get(i,1,k) + 2.0*gradKE.get(i,1,k); // + others //

        sl_state.aux1d->set(i,2,1, k, E     );
        sl_state.aux1d->set(i,2,2, k, Et    );
        sl_state.aux1d->set(i,2,3, k, Ett   );

    }

    ///////////////////////////////////////////////////////////////////////////
    // Everything below here is used for computing and saving Ettt:
    ///////////////////////////////////////////////////////////////////////////

    // terms used for 4th order stuff ... //
    // Ettt = -M3_xx + (2*E_x+rho)*M1 + 3*E*(M1)_x
    dTensorBC3 tmp0(mx,2,kmax1d,mbc,1);
    dTensorBC3 tmp1(mx,meqn,kmax1d,mbc,1);

    // compute the third moment and its 2nd derivative 
    dTensorBC3 M3    (mx,meqn,kmax1d,mbc,1);
    dTensorBC3 M3_x  (mx,2*meqn,kmax1d,mbc,1);
    IntegrateQ1dMoment3(q, M3 );
    SetPeriodicBndValues1D( M3 );
    FiniteDiff( 1, 2, M3, M3_x );

    // Compute 2*Ex+rho //
    FiniteDiff(1, 2*meqn, Enow, tmp0);
    for( int i=1; i <= mx; i++ )
    for( int k=1; k <= kmax1d; k++ )
    { tmp1.set(i,1,k, 2.0*tmp0.get(i,1,k) + rho.get(i,1,k) ); }

    // compute (2Ex+rho) * M1
    dTensorBC3 prod1(mx,meqn,kmax1d,mbc,1);
    MultiplyFunctions( tmp1, M1, prod1 );

    // compute M1_x and M1_xx
    dTensorBC3 M1_x  (mx,2*meqn,kmax1d,mbc,1);
    FiniteDiff( 1, 2, M1, M1_x );

    //compute 3*M1_x and E * (3*M1)_x
    dTensorBC3 prod2(mx,meqn,kmax1d,mbc,1);
    for( int i=1; i <= mx; i++ )
    for( int k=1; k <= kmax1d; k++ )
    { tmp1.set(i, 1, k, 3.0*M1_x.get(i, 1, k) ); }
    MultiplyFunctions( Enow, tmp1, prod2 );

    /////////////////////// Save Ettt /////////////////////////////////////////
    for( int i = 1; i <= mx; i++ )
    for( int k = 1; k <= kmax1d; k ++ )
    {
        double Ettt  = -M3_x.get(i,2,k) + prod1.get(i,1,k) + prod2.get(i,1,k);
        sl_state.aux1d->set(i,2,4, k, Ettt );
    }

}
///////////////////////////////////////////////////////////////////////////////
