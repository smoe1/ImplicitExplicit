#include <stdio.h>
#include <stdlib.h>
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"
//////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation forward in time.
//
// The advection equation: q_t + u q_x = 0 is solved by simply applying
// enough quadrature points at the future time, and tracing them backwards.
//
// This routine is the Non-Conservative "NonCons" form.
//
// It is assumed that aux(1,1,1) = u is constant.
//
//////////////////////////////////////////////////////////////////////////////
void StepAdvecNonCons(const double& dt, const dTensor2& node,
  dTensorBC3& auxvals, dTensorBC1& smax,
  const dTensorBC3& qold, dTensorBC3& qnew)
{

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

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points for the unit interval
    //////////////////////////////////////////////////////////////////////////

//  const int numpts = kmax + 2;
    const int numpts = 7;
    const double s_area = 2.0;
    dTensor1 wgt(numpts), spts(numpts);
    void setGaussLobattoPoints1d(dTensor1& x1d, dTensor1& wgt);
    setGaussLobattoPoints1d(spts, wgt);

    // compute quadrature points for where Q needs to be sampled
    const double eta = speed*dt/dx;
    dTensor1 spts_old ( numpts );
    iTensor1 ishift   ( numpts );
    void translateXi( 
        double eta, const dTensor1& xi, dTensor1 &xinew, iTensor1& ishift );
    translateXi( eta, spts, spts_old, ishift );

    // Legendre Basis functions, their derivatives, evaluated at all the
    // necessary poitns
    dTensor2 phi      (numpts, 5), phi_x      (numpts,5);
    dTensor2 phi_old  (numpts, 5), phi_x_old  (numpts,5);
    void evaluateLegendrePolys( 
            const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x);
    evaluateLegendrePolys( dx, spts, phi, phi_x);
    evaluateLegendrePolys( dx, spts_old, phi_old, phi_x_old);

#pragma omp parallel for
    for(int i=1; i<=melems; i++)
    {

        //loop over each equation
        for(int me=1; me<= meqn; me++)
        {

            for( int k=1; k <= kmax; k++)
            {
                double sum = 0.0;
                for( int m=1; m <= numpts; m++ )
                {

                    // periodicity is enforced here!
                    int io = (int)( i + ishift.get(m) );
                    io = iMod((io-1), melems)+1;

                    // evaluate q at the old point:
                    double qval = 0.0;
                    for( int kq=1; kq <= kmax; kq++ )
                    {
                        qval += qold.get(io, me, kq) * phi_old.get(m, kq);
                    }
                    sum += wgt.get(m) * qval * phi.get(m, k);
                }
                qnew.set(i, me, k, sum / s_area );
            }

        }//end of loop over each cell
    }//end of loop over each eqn  

    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    SetBndValues(node, auxvals, qnew);

}

///////////////////////////////////////////////////////////////////////////////
// This function sucks and shouldn't be used.  Although it works, ... it still
// sucks.
//
// This is being keps here for easy access to quadrature weights and points.
///////////////////////////////////////////////////////////////////////////////
/*
double GaussQuad(double a,double b, int numwghts, double (*Func)(double) )
{
    //local variables
    dTensor1 wgt(numwghts), spts(numwghts);
    double I = 0; //the value of the integral
    int i;

    //some error checking
    if( numwghts > 5 || numwghts < 1 ) {
        cout << "error in GaussQuad. ";
        cout << "bad integral inputs:" << endl;
        cout << "\t numwghts = " << numwghts << endl;
        cout << "\t a = " << a << '\t' << "b = " << b << endl;
        cout << "returning 0 " << endl;
        return 0;
    } else if( abs(a-b) < 1e-15) {
        return 0;
    } 

    // ---------------------------------
    // Set quadrature weights and points for the unit interval
    // ---------------------------------
    switch ( numwghts )
    {
        case 1:
          wgt.set(1,  2.0e0 );
      
          spts.set(1, 0.0e0 );
                        
          break;

        case 2:
          wgt.set(1,   1.0 );
          wgt.set(2,   1.0 );
      
          spts.set(1, -1.0/sq3 );
          spts.set(2,  1.0/sq3 );
                        
          break;

            case 3:
          wgt.set(1, 5.0e0/9.0e0 );
          wgt.set(2, 8.0e0/9.0e0 );
          wgt.set(3, 5.0e0/9.0e0 );
          
          spts.set(1, -sq3/sq5 );
          spts.set(2,  0.0e0 );
          spts.set(3,  sq3/sq5 );
                        
          break;

            case 4:
          wgt.set(1, (18.0 - sqrt(30.0))/36.0 );
          wgt.set(2, (18.0 + sqrt(30.0))/36.0 );
          wgt.set(3, (18.0 + sqrt(30.0))/36.0 );
          wgt.set(4, (18.0 - sqrt(30.0))/36.0 );
          
          spts.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
          spts.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
          spts.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
          spts.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );
                        
          break;

            case 5:
          wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
          wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
          wgt.set(3,  128.0/225.0 );
          wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
          wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );
          
          spts.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
          spts.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
          spts.set(3,  0.0 );
          spts.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
          spts.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
                        
          break;
    }

    //compute the integral sum approximation
    I = 0;    
    for( i=1; i <= numwghts; i++ )
    {
        I += wgt.get(i) * Func( (b-a)/2.0 * spts.get(i) + (a+b)/2.0 ) ;
    }
    return (b-a)/2.0 * I;

}
*/
