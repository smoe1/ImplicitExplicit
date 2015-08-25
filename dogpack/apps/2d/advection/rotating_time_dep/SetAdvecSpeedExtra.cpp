///////////////////////////////////////////////////////////////////////////////
// 
// 2d advection -- used to attain 4th order accuracy on
//
// q_t + u(y,t) q_x + v(x,t) q_t = 0;
// 
// By adding in extra terms into velocities u and v.
//
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "constants.h"
#include "SLState.h"


void SetAdvecSpeedExtra(
    const dTensorBC4& q, const dTensorBC4& aux, dTensorBC3& aux_extra,
    SL_state& sl_state )
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux_extra.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    const int kmax1d  = dogParams.get_space_order();
    const int mpoints = kmax1d*kmax1d;

    const double dtnow  = sl_state.dt_stages->get(sl_state.stage_num);
    const double tn     = sl_state.tn;

    const double t      = sl_state.t;
    const double dt2    = (t - tn) + dtnow;
    const double dt1    = (t - tn);

    // Mapping from 2D weights onto 1D weights
    iTensor1 k2d(5);

    // Here, u1 = u1(y):
    k2d.set(5,  15);
    k2d.set(4,  10);
    k2d.set(3,  6);
    k2d.set(2,  3);
    k2d.set(1,  1);

    // add in the weights to the appropriate spots
    for( int j = 1; j <= my; j++ )
    {
        for( int k = 1; k <= kmax1d; k ++ )
        {
            
            double E      = sl_state.aux1d->get(j,1,1,k);
            double Et     = sl_state.aux1d->get(j,1,2,k);
            double Ett    = sl_state.aux1d->get(j,1,3,k);
            double Ettt   = sl_state.aux1d->get(j,1,4,k);

            double  
                tmp  = E    * ( dtnow );
                tmp += Et   * ( pow(dt2,2) - pow(dt1,2)  ) / ( 2.0 );
                tmp += Ett  * ( pow(dt2,3) - pow(dt1,3)  ) / ( 6.0 );
                tmp += Ettt * ( pow(dt2,4) - pow(dt1,4)  ) / ( 24.0 );

            aux_extra.set(j, 1, k2d.get(k), tmp / dtnow );
        }
    }

    // Here, u2 = u2(x):
    k2d.set(5, 14);
    k2d.set(4,  9);
    k2d.set(3,  5);
    k2d.set(2,  2);
    k2d.set(1,  1);

    // add in the weights to the appropriate spots
    for( int i = 1; i <= mx; i++ )
    {
        for( int k = 1; k <= kmax1d; k ++ )
        {
            
            double E      = sl_state.aux1d->get(i,2,1,k);
            double Et     = sl_state.aux1d->get(i,2,2,k);
            double Ett    = sl_state.aux1d->get(i,2,3,k);
            double Ettt   = sl_state.aux1d->get(i,2,4,k);

            double  
                tmp  = E    * ( dtnow );
                tmp += Et   * ( pow(dt2,2) - pow(dt1,2)  ) / ( 2.0 );
                tmp += Ett  * ( pow(dt2,3) - pow(dt1,3)  ) / ( 6.0 );
                tmp += Ettt * ( pow(dt2,4) - pow(dt1,4)  ) / ( 24.0 );

            aux_extra.set(i, 2, k2d.get(k), tmp / dtnow );
        }
    }

}// end of function SetAdvecSpeedExtra
///////////////////////////////////////////////////////////////////////////////
