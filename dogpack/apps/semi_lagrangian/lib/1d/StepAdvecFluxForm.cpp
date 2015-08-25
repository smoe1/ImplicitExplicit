#include <stdio.h>
#include <stdlib.h>
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"
#include "dg_cell.h"


// return the shift index for xi //
inline int get_shift_index(double xi) 
{ return (int) floor( 0.5 * (xi+1.0) ); }

// shift xi to the interval [-1,1] //
inline double shift_xi(double xi) 
{ return xi - 2.0 * ( (double) get_shift_index(xi) ); }

//////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation forward in time.
//
// The advection equation: q_t + u q_x = 0.  
//
// It is assumed that aux(1,1,1) = u is constant.
//
//////////////////////////////////////////////////////////////////////////////
void StepAdvecFluxForm(const double& dt, const dTensor2& node,
  dTensorBC3& auxvals, dTensorBC1& smax,
  const dTensorBC3& qold, dTensorBC3& qnew)
{

    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    SetBndValues(node, auxvals, qnew);

    //-local parameters -----------------------------------------------------//
    const int melems  = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int kmax    = qnew.getsize(3);     //number of polynomials
    const int mbc     = qnew.getmbc();       //number of ghost cells
    //-----------------------------------------------------------------------//

    //save the center of each grid cell
    const double dx = node.get(2,1)-node.get(1,1);
    const double speed = auxvals.get(1,1,1);

    for(int j=1-mbc; j<=melems+mbc; j++)
    { smax.set(j, Max(smax.get(j), fabs(speed) ) ); }

    // compute quadrature points for where Q needs to be sampled
    const double s_area = 2.0;
    const double large_eta = speed*dt/dx;

    // integration points and weights (for interior)
    const int mpts = 5;

    dTensor1 wgt( mpts ), spts( mpts );
    dTensor2 phi( mpts, 5), phi_x(mpts, 5);
    if( kmax > 1 )
    {
        void setGaussPoints1d(dTensor1& x1d, dTensor1& w1d);
        void evaluateLegendrePolys( 
                const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x);

        setGaussPoints1d( spts, wgt );
        evaluateLegendrePolys(dx, spts, phi, phi_x);
    }

#pragma omp parallel for
    for(int i=1; i<=melems; i++)
    {

        //loop over each equation
        for(int me=1; me<= meqn; me++)
        {

            dg_cell* dgc = new dg_cell();

            for( int k=1; k <= kmax; k++)
            {

                double qnp1 = qold.get(i,me,k);

                if( k > 1 )
                {
                    // interior integral
                    for( int m=1; m <= mpts; m++ )
                    {

                        double a   = spts.get(m) - 2.0 * large_eta;
                        double b   = spts.get(m);

                        int ishift = get_shift_index( a );

                        double integral = 0.0;
                        if( ishift == 0 )
                        {
                            // no shift takes place, can simply integrate from a to b
                            for( int kq = 1; kq <= kmax; kq++ )
                            {
                                integral += qold.get(i, me, kq) * 
                                    dgc->integratePhi( a, b, kq );
                            }

                        }
                        else if( ishift < 0 )
                        {

                            // integral at the shifted value
                            for( int kq = 1; kq <= kmax; kq++ )
                            {
                                integral += qold.get(i+ishift, me, kq) * 
                                    dgc->integratePhi( shift_xi(a), 1.0, kq );
                            }

                            // integrals between i and ishift
                            for( int n=-1; n > ishift; n-- )
                            { integral += qold.get(i + n, meqn, 1) * dgc->integratePhi(-1.0, 1.0, 1); }

                            // integral inside current cell
                            for( int kq = 1; kq <= kmax; kq++ )
                            {
                                integral += qold.get(i, me, kq) * 
                                    dgc->integratePhi( -1.0, spts.get(m), kq );
                            }

                        }
                        else
                        {
                            printf("   StepAdvecFluxForm.cpp problem can't handle negative velocities yet\n");
                            exit(1);
                        }

                        qnp1 += 0.25 * dx * wgt.get(m) * phi_x.get(m,k) * integral;

                    }
                }

                // left flux
                double a = -1.0 - 2.0 * large_eta;
                int ishift = get_shift_index( a );

                double Fm = 0.0;
                if( ishift < 0 )
                {
                    for( int kq = 1; kq <= kmax; kq++ )
                    {
                        Fm += qold.get(i+ishift, me, kq) * 
                            dgc->integratePhi( shift_xi(a), 1.0, kq ) / 2.0;
                    }
                    for( int n=-1; n > ishift; n-- )
                    {
                        Fm += qold.get(i+n, me, 1);
                    }
                }
                Fm *= sqrt(2*k-1) * pow(-1, k+1 );

                // right flux
                double Fp = 0.0;
                a += 2.0;
                ishift = get_shift_index( a );
                if( ishift < 1 )
                {
                    for( int kq = 1; kq <= kmax; kq++ )
                    {
                        Fp += qold.get(i+ishift, me, kq) * 
                            dgc->integratePhi( shift_xi(a), 1.0, kq ) / 2.0;
                    }
                }
                for( int n=0; n > ishift; n-- )
                { Fp += qold.get(i+n, me, 1); }
                Fp *= sqrt(2*k-1);

                qnew.set(i, me, k, qnp1 - ( Fp - Fm ) );

            }

            delete dgc;

        }//end of loop over each cell
    }//end of loop over each eqn  

    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    SetBndValues(node, auxvals, qnew);

}
