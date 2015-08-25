#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

///////////////////////////////////////////////////////////////////////////////
// Routine for setting advection speeds for advection equation of the form
// q_t + u(y) q_x + v(x) q_y = 0.
//
// Function also sets the maximum wave speed for each cell.
//
//  Note: for mpoints1d = 3, mpoints are arranged as follows:
//   
//   ________________
//   |              |
//   |  7    8   9  |
//   |              |
//   |  4    5   6  |
//   |              | 
//   |  1    2   3  |
//   |______________|
//
// Params:
//    aux(mx, my, maux, kmax, mbc) - maux = 1 is u(y), maux = 2 is v(x)
//
//    The advection speeds sampled at the gaussian quadrature points 
//    1:mpoints1d are provided in:
//
//    u1(my, mpoints1d) = u(y) - the 1d gauss point values for speeds
//    u2(mx, mpoints1d) = v(x) - the 1d gauss point values for speeds
//
///////////////////////////////////////////////////////////////////////////////
void SetAdvecSpeed(const dTensor2& phi, const dTensorBC4& q,
		   const dTensorBC4& aux, dTensorBC3& smax, const int &vel_dir,
		   dTensorBC2& u1, dTensorBC2& u2, SL_state& sl_state)
{

    const int mx   = aux.getsize(1);
    const int my   = aux.getsize(2);
    const int maux = aux.getsize(3);
    const int kmax = aux.getsize(4);
    const int mbc  = aux.getmbc();

    const int mpoints1d = u1.getsize(2);
    const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    const int ndims = 2;
    dTensorBC3 aux_extra( Max(mx,my), ndims, kmax, mbc );

    for( int j = 1; j <= my; j++ )
    {
        for( int k = 1; k <= kmax; k++ )
        {
            aux_extra.set(j, 1, k, aux.get(1,j,1,k) );
        }
    }

    for( int i = 1; i <= mx; i++ )
    {
        for( int k = 1; k <= kmax; k++ )
        {
            aux_extra.set(i, 2, k, aux.get(i,1,2,k) );
        }
    }

    // Applications that require extra coefficients are required to overrride
    // this:
    //
    // For the case of VP + Yoshida splitting, we only add in extra terms when
    // advecting in the direction of even stage numbers
    //
    // For this problem, we add in the corrections every time step:
    void SetAdvecSpeedExtra( const dTensorBC4& q, const dTensorBC4& aux, 
            dTensorBC3& aux_extra, SL_state& sl_state );
    { SetAdvecSpeedExtra(q, aux, aux_extra, sl_state ); }

    switch( vel_dir )
    {
        case 1:

            // save the velocity u1
            for(int j=1; j <= my; j++)
            {
                for(int m1d = 1; m1d<=mpoints1d; m1d++)
                {

                    // TODO -- this is wasted space!!!  there's no need to sum
                    // over all kmax polynomials, since over half of them are
                    // zero
                    double tmp1 = 0.0;
                    for(int k=1; k<=kmax; k++)
                    {
                        // sample first row 
                        tmp1 += aux_extra.get(j,1,k)*phi.get( (m1d-1)*mpoints1d+1, k);
                    }
                    u1.set(j, m1d, tmp1);                // velocity, u(y)
                }
            }
            break;

        case 2:

            // save the velocity u2
            //for(int i=(1-mbc); i <= mx+mbc; i++)
            for(int i=1; i <= mx; i++)
            {
                for(int m1d = 1; m1d<=mpoints1d; m1d++)
                {

                    // TODO -- this is wasted space!!!  there's no need to sum
                    // over all kmax polynomials, since over half of them are
                    // zero
                    double tmp2 = 0.0;
                    for(int k=1; k<=kmax; k++)
                    {
                        // sample first collumn:
                        tmp2 += aux_extra.get(i,2,k)*phi.get( m1d,  k);
                    }

                    u2.set(i,m1d,tmp2);                 // velocity, v(x)
                }

            }
            break;
    }

    // set maximum wave speed
    // TODO - just do this along one row?
    const int m = (int)((mpoints1d+1)/2);
    for(int i=(1); i<= mx; i++)
    for(int j=(1); j<= my; j++)
    {

        double tmp1 = 0.0;
        double tmp2 = 0.0;
        {
            tmp1 = Max(tmp1, fabs( dy*u1.get(j,m) ) );
            tmp2 = Max(tmp2, fabs( dx*u2.get(i,m) ) );    
        }
        smax.set(i,j,1, tmp1 );
        smax.set(i,j,2, tmp2 );
    }

}// end of function SetAdvecSpeed
///////////////////////////////////////////////////////////////////////////////
