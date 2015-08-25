///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

void SetPeriodicBndValues1D( dTensorBC3& q);
void ConvertQ2dToQ1d(const int &mopt, int istart, int iend, 
                     const dTensorBC4& qin, dTensorBC3& qout);

// The purpose of this routine is to compute the necessary terms for the
// expansion in 
//
//    \bar{E}(t) = E^n + (t-tn) E_t + (t-tn)^2 / 2 * E_tt + (t-tn)^3 / 3! E_ttt
//
// Results are saved into sl_state->aux1d.  See SLState.h for what that
// auxiliary array represents
void InitSLState( 
    const dTensorBC4& q, const dTensorBC4& aux, SL_state& sl_state )
{

    // We only ever use the expansion when running 4th order code
    if( dogParams.get_time_order() < 4 )
    { return; }

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    const int kmax1d  = dogParams.get_space_order();
    const int mpoints = kmax1d*kmax1d;

    const double tn     = sl_state.tn;

    // Mapping from 2D weights onto 1D weights
    iTensor1 k2d(5);

    // Here, u1 = u1(v):
    k2d.set(5,  15);
    k2d.set(4,  10);
    k2d.set(3,  6);
    k2d.set(2,  3);
    k2d.set(1,  1);

    // -- Velocity in the x-direction including time derivatives: -- //
    for( int j = 1; j <= my; j++ )
    for( int k1 = 1; k1 <= kmax1d; k1++ )
    {

        // k - derivatives of the velocity:
        double u    = exp(tn) * aux.get(1, j, 1, k2d.get(k1) );
        double ut   = exp(tn) * aux.get(1, j, 1, k2d.get(k1) );
        double utt  = exp(tn) * aux.get(1, j, 1, k2d.get(k1) );
        double uttt = exp(tn) * aux.get(1, j, 1, k2d.get(k1) );

        sl_state.aux1d->set(j, 1, 1, k1, u     );
        sl_state.aux1d->set(j, 1, 2, k1, ut   );
        sl_state.aux1d->set(j, 1, 3, k1, utt  );
        sl_state.aux1d->set(j, 1, 4, k1, uttt );

    }

    // Here, u2 = u2(x):
    k2d.set(5, 14);
    k2d.set(4,  9);
    k2d.set(3,  5);
    k2d.set(2,  2);
    k2d.set(1,  1);

    // -- Velocity in the y-direction including time derivatives: -- //
    for( int i = 1; i <= mx; i++ )
    for( int k1 = 1; k1 <= kmax1d; k1++ )
    {

        // k - derivatives of the velocity:
        double u    = exp(tn) * aux.get(i,1, 2, k2d.get(k1) );
        double ut   = exp(tn) * aux.get(i,1, 2, k2d.get(k1) );
        double utt  = exp(tn) * aux.get(i,1, 2, k2d.get(k1) );
        double uttt = exp(tn) * aux.get(i,1, 2, k2d.get(k1) );

        sl_state.aux1d->set(i, 2, 1, k1, u     );
        sl_state.aux1d->set(i, 2, 2, k1, ut   );
        sl_state.aux1d->set(i, 2, 3, k1, utt  );
        sl_state.aux1d->set(i, 2, 4, k1, uttt );

    }

}
///////////////////////////////////////////////////////////////////////////////
