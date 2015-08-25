//////////////////////////////////////////////////////////////////////////////
// Function to add in correction terms to aux_extra
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "constants.h"
#include "SLState.h"

// derivative of electric field, E, E_t
//void ElectricField(const dTensor2& xpts, const dTensor2& qvals, 
//        const dTensor2& auxvals, dTensor2& e, void* data)
void ElectricField(const dTensor2& xpts, dTensor2& e, void* data)
{
    const int numpts=xpts.getsize(1);

    // grab the current time and the next time
    double t = *(double*)data;

    double E, Et, Ett, Ettt;

    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i,1);
        double v = xpts.get(i,2);

        E    = -sqrt(pi) / 4.0 * sin( 2.0*(x-pi*t) );
        Et   = -sqrt(pi) / 4.0 * cos( 2.0*(x-pi*t) ) * (-2.0*pi);
        Ett  =  sqrt(pi) / 4.0 * sin( 2.0*(x-pi*t) ) * pow((-2.0*pi),2);
        Ettt =  sqrt(pi) / 4.0 * cos( 2.0*(x-pi*t) ) * pow((-2.0*pi),3);

        e.set(i, 1, E    );
        e.set(i, 2, Et   );
        e.set(i, 3, Ett  );
        e.set(i, 4, Ettt );

    }

}

// extra source terms involved in E, Et, Ett, Ettt //
//void ExtraSource(const dTensor2& xpts, const dTensor2& qvals, 
//        const dTensor2& auxvals, dTensor2& e, void* data)
void ExtraSource(const dTensor2& xpts, 
        dTensor2& e, void* data)
{
    const int numpts=xpts.getsize(1);
    double x, v, t;

    // grab the current time and the next time
    t = *(double*)data;

    for (int i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        v = xpts.get(i,2);

        // these are direct functions of the source term ONLY //

        // Et = -rho_u + g //
        double g = (1.0/8.0)*sqrt(pi)*(4.0*cos(-2.0*x+2.0*pi*t)*pi+2.0-cos(-2.0*x+2.0*pi*t));

        // Ett = -rho * E + \partial_x ( 2*KE ) + g1 //
        double g1 =
            -(1.0/16.0)*sin(-2*x+2*pi*t)*(16.0*pow( pi,(5.0/2.0) ) -
            4.0*pi+2.0*cos(-2.0*x+2.0*pi*t)*pi-3.0*sqrt(pi));

        // Ettt = (rho-2)*rho_u + E (rho_u)_x + M3_x + g2 //
        //         where M3 = \int v^3 f dv
        double g2 =
        -2.0*pow( pi, (7.0/2.0) )*cos(-2.0*x+2.0*pi*t)+(7.0/32.0)*sqrt(pi)*cos(-2.0*x+2.0*pi*t)+(1.0/2.0)*cos(-2.0*x+2.0*pi*t)*pi-(3.0/8)*pi*pow(
        cos(-2.0*x+2.0*pi*t),2.0 ) -(1.0/16.0)*pi;

        e.set(i, 1, 0.0  );         // added to E.
        e.set(i, 2, g    );         // added to -rho_u.  Et  = -rho_u + g(x,t)
        e.set(i, 3, g1 );           // added to Ett  = -2*KE_x - \rho * E + g1(x,t).     
        e.set(i, 4, g2 );           // added to Ettt=

    }

}

void SetAdvecSpeedExtra( 
    const dTensorBC4& q, const dTensorBC4& aux,
    dTensorBC3& aux_extra, SL_state& sl_state )
{

    const int space_order = dogParams.get_space_order();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux_extra.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    const int mpoints = space_order*space_order;
    const int kmax1d  = space_order;

    const double dtnow  = sl_state.dt_stages->get(sl_state.stage_num);
    const double tn     = sl_state.tn;

    const double t      = sl_state.t;
    const double dt2    = (t - tn) + dtnow;
    const double dt1    = (t - tn);

    iTensor1 k2d(5);
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
