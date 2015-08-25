#include <cmath>
#include <iostream>
#include "dogdefs.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//     
//     Function to solve the Poisson equation: 
//
//                u_xx = -f(x) on [a,b]
//                u'(a) = gamma    u(b) = beta
//
//     This equation is solved via a discontinuous Galerkin (DG) method
//     of order 1,2,3,4, or 5. In the DG method, we actually solve the
//     Poisson equation by writing it as a system of two first order
//     differential equations:
//
//                v_x = f(x) on [a,b]
//                u_x = -v   on [a,b]
//                v(a) = -gamma    u(b) = beta
//
//     Combinations of other boundary conditions u(a) = alpha; u'(b) = delta
//     Can be specified by first calling /PoissonSolver/ConvertBC.cpp
//
//     Parameters:
//
//     Fvals(melems, 1, morder, mbc) - Legendre weights for right-hand     
//                                     side function f(x).  Note: this carries
//                                     a minus sign with it!
//                    
//     qvals(melems, meqn, morder, mbc) - legendre weights for solution u.
//         qvals(:,mstart,:)   = 'electric field'  .... the gradient, -u_x = v
//         qvals(:,mstart+1,:) = 'potential' ... the solution u(x)
//
///////////////////////////////////////////////////////////////////////////////

void PoissonSolve(const int mstart, const dTensor2& node, const double gamma, 
        const double beta, const dTensorBC3& Fvals, dTensorBC3& qvals)
{
    //parameters
    int i,k,n,l;
    int melems = qvals.getsize(1);
    int meqn   = qvals.getsize(2);
    int morder = qvals.getsize(3);
    double count,tmp;
    double dx = node.get(2,1)-node.get(1,1);

    //matrices and inverse matrices needed are Ainv and Cinv
    dTensor2 Ainv(morder,morder);
    dTensor2 Cinv(morder,morder);
    dTensor2 Bm(5,5),Dm(5,5);
    dTensor1 tmpVec(morder);    //helper for doing matrix multiplication

    // quick error check
    if ( (mstart+1)>meqn )
    {
        cout << " ERROR in PoissonSolve: mstart is incorrectly set" << endl;
        cout << "    mstart = " << mstart << endl;
        cout << "      meqn = " << meqn << endl;
        exit(1);
    }


    // ------------------------------------------------------------------
    //
    //     SETUP BLOCK MATRICES, A^{-1}, C^{-1}, B, AND D
    //
    // ------------------------------------------------------------------    
    switch(morder)
    {      
        case 1:
            Ainv.set(1,1, 1.0 );
            Cinv.set(1,1, 1.0 );
            break;

        case 2:
            Ainv.set(1,1,  0.5 );
            Ainv.set(1,2, -sq3/6.0 );
            Ainv.set(2,1,  sq3/6.0 );
            Ainv.set(2,2,  1.0/6.0 );

            Cinv.set(1,1,  0.5 );
            Cinv.set(1,2,  sq3/6.0 );
            Cinv.set(2,1, -sq3/6.0 );
            Cinv.set(2,2,  1.0/6.0 );
            break;

        case 3:	
            Ainv.set(1,1,  0.5 );
            Ainv.set(1,2, -sq3/6.0 );
            Ainv.set(1,3,  0.0 );
            Ainv.set(2,1,  sq3/6.0 );
            Ainv.set(2,2,  0.0 );
            Ainv.set(2,3, -sq3*sq5/30.0 );
            Ainv.set(3,1,  0.0 );
            Ainv.set(3,2,  sq3*sq5/30.0 );
            Ainv.set(3,3,  0.1 );

            Cinv.set(1,1,  0.5 );
            Cinv.set(1,2,  sq3/6.0 );
            Cinv.set(1,3,  0.0 );
            Cinv.set(2,1, -sq3/6.0 );
            Cinv.set(2,2,  0.0 );
            Cinv.set(2,3,  sq3*sq5/30.0 );
            Cinv.set(3,1,  0.0 );
            Cinv.set(3,2, -sq3*sq5/30.0 );
            Cinv.set(3,3,  0.1 );	
            break;

        case 4:	  

            Ainv.set(1,1,  0.5 );
            Ainv.set(1,2, -sq3/6.0 );
            Ainv.set(1,3,  0.0 );
            Ainv.set(1,4,  0.0 );
            Ainv.set(2,1,  sq3/6.0 );
            Ainv.set(2,2,  0.0 );
            Ainv.set(2,3, -sq3*sq5/30.0 );
            Ainv.set(2,4,  0.0 );
            Ainv.set(3,1,  0.0 );
            Ainv.set(3,2,  sq3*sq5/30.0 );
            Ainv.set(3,3,  0.0 );
            Ainv.set(3,4, -sq5*sq7/70.0 );
            Ainv.set(4,1,  0.0 );
            Ainv.set(4,2,  0.0 );
            Ainv.set(4,3,  sq5*sq7/70.0 );
            Ainv.set(4,4,  1.0/14.0 );

            Cinv.set(1,1,  0.5 );
            Cinv.set(1,2,  sq3/6.0 );
            Cinv.set(1,3,  0.0 );
            Cinv.set(1,4,  0.0 );
            Cinv.set(2,1, -sq3/6.0 );
            Cinv.set(2,2,  0.0 );
            Cinv.set(2,3,  sq3*sq5/30.0 );
            Cinv.set(2,4,  0.0 );
            Cinv.set(3,1,  0.0 );
            Cinv.set(3,2, -sq3*sq5/30.0 );
            Cinv.set(3,3,  0.0 );
            Cinv.set(3,4,  sq5*sq7/70.0 );
            Cinv.set(4,1,  0.0 );
            Cinv.set(4,2,  0.0 );
            Cinv.set(4,3, -sq5*sq7/70.0 );
            Cinv.set(4,4,  1.0/14.0 );	
            break;

        case 5:

            Ainv.set(1,1,  0.5 );
            Ainv.set(1,2, -sq3/6.0 );
            Ainv.set(1,3,  0.0 );
            Ainv.set(1,4,  0.0 );
            Ainv.set(1,5,  0.0 );
            Ainv.set(2,1,  sq3/6.0 );
            Ainv.set(2,2,  0.0 );
            Ainv.set(2,3, -sq3*sq5/30.0 );
            Ainv.set(2,4,  0.0 );
            Ainv.set(2,5,  0.0 );
            Ainv.set(3,1,  0.0 );
            Ainv.set(3,2,  sq3*sq5/30.0 );
            Ainv.set(3,3,  0.0 );
            Ainv.set(3,4, -sq5*sq7/70.0 );
            Ainv.set(3,5,  0.0 );
            Ainv.set(4,1,  0.0 );
            Ainv.set(4,2,  0.0 ); 
            Ainv.set(4,3,  sq5*sq7/70.0 );
            Ainv.set(4,4,  0.0 );
            Ainv.set(4,5, -sq7/42.0 );
            Ainv.set(5,1,  0.0 );
            Ainv.set(5,2,  0.0 );
            Ainv.set(5,3,  0.0  );
            Ainv.set(5,4,  sq7/42.0 );
            Ainv.set(5,5,  1.0/18.0 );

            Cinv.set(1,1,  0.5 );
            Cinv.set(1,2,  sq3/6.0 );
            Cinv.set(1,3,  0.0 );
            Cinv.set(1,4,  0.0 );
            Cinv.set(1,5,  0.0 );
            Cinv.set(2,1, -sq3/6.0 );
            Cinv.set(2,2,  0.0 );
            Cinv.set(2,3,  sq3*sq5/30.0 );
            Cinv.set(2,4,  0.0 );
            Cinv.set(2,5,  0.0 );
            Cinv.set(3,1,  0.0);
            Cinv.set(3,2, -sq3*sq5/30.0 );
            Cinv.set(3,3,  0.0);
            Cinv.set(3,4,  sq5*sq7/70.0 );
            Cinv.set(3,5,  0.0 );
            Cinv.set(4,1,  0.0 );
            Cinv.set(4,2,  0.0 );
            Cinv.set(4,3, -sq5*sq7/70.0 );
            Cinv.set(4,4,  0.0 );
            Cinv.set(4,5,  sq7/42.0 );
            Cinv.set(5,1,  0.0 );
            Cinv.set(5,2,  0.0 );
            Cinv.set(5,3,  0.0 );
            Cinv.set(5,4, -sq7/42.0 );
            Cinv.set(5,5,  1.0/18.0 );
            break;
    }//end of writing matrices Ainv, Cinv

    Bm.set(1,1, -1.0 );
    Bm.set(1,2, -sq3 );
    Bm.set(1,3, -sq5 );
    Bm.set(1,4, -sq7 );
    Bm.set(1,5, -3.0);
    Bm.set(2,1,  sq3 );
    Bm.set(2,2,  3.0);
    Bm.set(2,3,  sq3*sq5 );
    Bm.set(2,4,  sq3*sq7 );
    Bm.set(2,5,  3.0*sq3 );
    Bm.set(3,1, -sq5 );
    Bm.set(3,2, -sq3*sq5 );
    Bm.set(3,3, -5.0 );
    Bm.set(3,4, -sq5*sq7 );
    Bm.set(3,5, -3.0*sq5 );
    Bm.set(4,1,  sq7 );
    Bm.set(4,2,  sq3*sq7 );
    Bm.set(4,3,  sq5*sq7 );
    Bm.set(4,4,  7.0 );
    Bm.set(4,5,  3.0*sq7 );
    Bm.set(5,1, -3.0 );
    Bm.set(5,2, -3.0*sq3 );
    Bm.set(5,3, -3.0*sq5 );
    Bm.set(5,4, -3.0*sq7 );
    Bm.set(5,5, -9.0 );

    Dm.set(1,1, -1.0 );
    Dm.set(1,2,  sq3 );
    Dm.set(1,3, -sq5 );
    Dm.set(1,4,  sq7 );
    Dm.set(1,5, -3.0 );
    Dm.set(2,1, -sq3 );
    Dm.set(2,2,  3.0);
    Dm.set(2,3, -sq3*sq5 );
    Dm.set(2,4,  sq3*sq7 );
    Dm.set(2,5, -3.0*sq3 );
    Dm.set(3,1, -sq5 );
    Dm.set(3,2,  sq3*sq5 );
    Dm.set(3,3, -5.0 );
    Dm.set(3,4,  sq5*sq7 );
    Dm.set(3,5, -3.0*sq5 );
    Dm.set(4,1, -sq7 );
    Dm.set(4,2,  sq3*sq7 );
    Dm.set(4,3, -sq5*sq7 );
    Dm.set(4,4,  7.0 );
    Dm.set(4,5, -3.0*sq7 );
    Dm.set(5,1, -3.0 );
    Dm.set(5,2,  3.0*sq3 );
    Dm.set(5,3, -3.0*sq5 );
    Dm.set(5,4,  3.0*sq7 );
    Dm.set(5,5, -9.0 );


    // ------------------------------------------------------------------
    //
    //     FORWARD SUBSTITUTION TO SOLVE FOR  ``ELECTRIC FIELD''
    //
    // ------------------------------------------------------------------

    //special case for the first vector, qvals(1,mstart,:)
    for(k=1; k<=morder; k++) 
    {
        if( k%2 == 0 )
        {
            tmpVec.set(k, dx*Fvals.get(1,1,k) + sqrt(2.0*k-1.0)*gamma);
        }
        else
        {
            tmpVec.set(k, dx*Fvals.get(1,1,k) - sqrt(2.0*k-1.0)*gamma);
        }
    }

    // qvals(1,mstart,:) = invA * tmpVec
    for(k=1; k<=morder; k++)
    {
        tmp = 0.0;
        for(n=1; n<=morder; n++)
        {  tmp = tmp + Ainv.get(k,n)*tmpVec.get(n);  }

        qvals.set(1,mstart,k, tmp );

    }

    //every other case, qvals(i,mstart,:) 
    //                   = invA*( dx*Fvals(i,1,:) - B*qvals(i-1,mstart,:) )
    for(i=2; i<=melems; i++)
    {
        //compute B*qvals(i-1,mstart,:)
        //save as tmpVec
        for(l=1; l<=morder; l++ )
        {
            tmp = 0.0;
            for(k=1; k<=morder; k++)
            {  tmp = tmp + Bm.get( l , k ) * qvals.get(i-1,mstart,k);  }

            tmpVec.set( l , tmp );
        }

        //compute dx*Fvals(i,1,:) - tmpVec
        //save result as tmpVec
        for(n=1; n<=morder; n++)
        {  tmpVec.set( n , dx*Fvals.get( i,1, n ) - tmpVec.get(n) );  }

        //compute Ainv * ( dx*Fvals(i,1,:) - tmpVec ) = Ainv * tmpVec
        for(l=1; l<=morder; l++)
        {
            tmp = 0;
            for(n=1; n<=morder; n++)
            {  tmp = tmp + Ainv.get(l, n ) * tmpVec.get(n);  }
            qvals.set(i,mstart,l, tmp );
        }

    }


    // ------------------------------------------------------------------
    //
    //     BACKWARD SUBSTITUTION TO SOLVE FOR  ``POTENTIAL''
    //
    // ------------------------------------------------------------------

    //special case for first substitution, qvals(melems,mstart+1,:)
    for(l=1; l<=morder; l++)
    {  tmpVec.set(l, sqrt(2.0*l-1.0)*beta + dx*qvals.get(melems,mstart,l) );  }

    for(k=1; k<=morder; k++)
    {
        count = 0.0;
        for(n=1; n<=morder; n++)
        {  count = count + Cinv.get(k,n)*tmpVec.get(n);  }
        qvals.set(melems,mstart+1,k, count);
    }

    //every other case, qvals(i,mstart+1,:) 
    //          = invC*( dx*qvals(i,mstart,:) - C*qvals(i+1,mstart+1,:) )
    for(i=(melems-1); i>0; i--)
    {
        //compute D*qvals(i+1,mstart+1,:)
        //save as tmpVec
        for(n=1; n<=morder ; n++)
        {
            tmp = 0.0;
            for(k=1; k<=morder; k++)
            {  tmp = tmp + Dm.get( n , k ) * qvals.get(i+1,mstart+1,k);  }
            tmpVec.set( n , tmp );
        }

        //compute dx*qvals(i,mstart,:) - tmpVec
        //save result as tmpVec
        for(n=1; n<=morder; n++)
        {  tmpVec.set( n , dx*qvals.get(i,mstart,n) - tmpVec.get(n) );  }

        //compute Cinv * ( dx*qvals(i,mstart,:) - tmpVec ) = Cinv * tmpVec
        for(l=1; l<=morder; l++)
        {
            tmp = 0.0;
            for(n=1; n<=morder; n++)
            {  tmp = tmp + Cinv.get(l, n ) * tmpVec.get(n);  }
            qvals.set(i,mstart+1,l, tmp );
        }
    }

}// end of function PoissonSolve.cpp
