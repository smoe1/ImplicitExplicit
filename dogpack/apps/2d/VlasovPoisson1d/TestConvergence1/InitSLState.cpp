///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

void FiniteDiff(int mestart, int meend, const dTensorBC3& q, dTensorBC3& qout);
void SetBndValues1D( dTensorBC3& q);
void L2Project_extra(const int istart,
                     const int iend,
                     const int jstart,
                     const int jend,
                     const int QuadOrder,
                     const int BasisOrder_qin,
                     const int BasisOrder_auxin,
                     const int BasisOrder_fout,    
                     const dTensorBC4* qin,
                     const dTensorBC4* auxin,
                     dTensorBC4* fout,
                     void (*Func)(const dTensor2&,const dTensor2&,
                                  const dTensor2&,dTensor2&, void* data),
                     void* data);

void ExtraSourceWrap(const dTensor2& xpts, 
		     const dTensor2& NOT_USED_1, 
		     const dTensor2& NOT_USED_2,
		     dTensor2& e, 
		     void* data)
{
  void ExtraSource(const dTensor2& xpts, 		   
		   dTensor2& e, 
		   void* data);

  ExtraSource(xpts,e,data);
}

void InitSLState( 
    const dTensorBC4& q, const dTensorBC4& aux, SL_state& sl_state )
{

    void ComputeElecField(double t, const dTensor2& node1d,
            const dTensorBC4& qvals, dTensorBC3& Evals);
    void ConvertQ2dToQ1d(const int &mopt, int istart, int iend, 
                     const dTensorBC4& qin, dTensorBC3& qout);
    void IntegrateQ1dMoment1(const dTensorBC4& q2d, dTensorBC3& q1d);

    const int space_order = dogParams.get_space_order();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int mpoints = space_order*space_order;
    const int kmax1d  = space_order;

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
    void IntegrateQ1dMoment2(const dTensorBC4& q2d, dTensorBC3& q1d);
    IntegrateQ1dMoment2(q, KE);  // 1/2 * \int v^2 * f dv //

    SetBndValues1D( KE );
    SetBndValues1D( Enow );
    SetBndValues1D( M1 );

    // do something to compute KE_x ... //
    dTensorBC3 gradKE(mx,2,kmax1d,mbc,1);
    FiniteDiff( 1, 2, KE, gradKE ); // compute KE_x and KE_xx //

    dTensorBC3 rho(mx,meqn,kmax1d,mbc,1);
    dTensorBC3 prod(mx,meqn,kmax1d,mbc,1);
    void IntegrateQ1d(const int mopt, const dTensorBC4& q2d, dTensorBC3& q1d);
    void MultiplyFunctions(const dTensorBC3& Q1, const dTensorBC3& Q2,
        dTensorBC3& qnew);
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

    // terms used for 4th order stuff ... //
    // Ettt = -M3_xx + (2*E_x+rho)*M1 + 3*E*(M1)_x
    dTensorBC3 tmp0(mx,2*meqn,kmax1d,mbc,1);
    dTensorBC3 tmp1(mx,meqn,kmax1d,mbc,1);

    // compute the third moment and its 2nd derivative 
    dTensorBC3 M3    (mx,meqn,kmax1d,mbc,1);
    dTensorBC3 M3_x  (mx,2*meqn,kmax1d,mbc,1);
    void IntegrateQ1dMoment3(const dTensorBC4& q2d, dTensorBC3& q1d);
    IntegrateQ1dMoment3(q, M3 );
    SetBndValues1D( M3 );
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
    ///////////////////////////////////////////////////////////////////////////

    if ( dogParams.get_source_term()>0 )
    {


        double* t_ptr = new double;
        *t_ptr = tn;

        // 2D Electric field, and 1D Electric Field
        dTensorBC4* ExactE;
        ExactE     = new dTensorBC4(mx,1,4,kmax,mbc);

        dTensorBC3* ExactE1d;
        ExactE1d   = new dTensorBC3(mx,4,kmax1d,mbc,1);

        // Extra Source terms needed for the electric field
        dTensorBC4* ExtraE;
        ExtraE     = new dTensorBC4(mx,1,4,kmax,mbc);

        dTensorBC3* ExtraE1d;
        ExtraE1d   = new dTensorBC3(mx,4,kmax1d,mbc,1);


        // Exact, electric field: (TODO - REMOVE THIS!)
//      void ElectricField(const dTensor2& xpts, dTensor2& e, void* data);
//      L2Project_extra(1-mbc, mx+mbc, 1, 1, space_order, -1, ExactE, 
//                      &ElectricField, (void*)t_ptr );
//      ConvertQ2dToQ1d(1, 1, mx, *ExactE, *ExactE1d);

        // Extra Source terms:
        void ExtraSourceWrap(const dTensor2& xpts, 
			     const dTensor2& NOT_USED_1, 
			     const dTensor2& NOT_USED_2,
			     dTensor2& e, 
			     void* data);
        L2Project_extra(1, mx, 1, 1, 20, space_order,
			space_order, space_order,
			&q, &aux, ExtraE, &ExtraSourceWrap, (void*)t_ptr );
        ConvertQ2dToQ1d(1, 1, mx, *ExtraE, *ExtraE1d);
  
        for( int i=1; i <= mx; i++ )
        for( int k=1; k <= kmax1d; k++ )
        {

            // electric fields w/o source term parts added in //
            double Et    = sl_state.aux1d->get(i,2,2,k);
            double Ett   = sl_state.aux1d->get(i,2,3,k);
            double Ettt  = sl_state.aux1d->get(i,2,4,k); 

            // add in missing terms from previously set values //
            sl_state.aux1d->set(i,2,2, k, Et    + ExtraE1d->get(i,2,k) );
            sl_state.aux1d->set(i,2,3, k, Ett   + ExtraE1d->get(i,3,k) );
            sl_state.aux1d->set(i,2,4, k, Ettt  + ExtraE1d->get(i,4,k) );

//          sl_state.aux1d->set(i,2,1, k, 0. );
//          sl_state.aux1d->set(i,2,2, k, 0. );
//          sl_state.aux1d->set(i,2,3, k, 0. );
//          sl_state.aux1d->set(i,2,4, k, 0. );


            // ADD IN EXACT ELECTRIC FIELD HERE:
//          sl_state.aux1d->set(i,2,1, k, ExactE1d->get(i,1,k) );
//          sl_state.aux1d->set(i,2,2, k, ExactE1d->get(i,2,k) );
//          sl_state.aux1d->set(i,2,3, k, ExactE1d->get(i,3,k) );
//          sl_state.aux1d->set(i,2,4, k, ExactE1d->get(i,4,k) );


        }

        delete ExactE;
        delete ExactE1d;
        delete t_ptr;
        delete ExtraE;
        delete ExtraE1d;

    }

    
}
///////////////////////////////////////////////////////////////////////////////

// Periodic boundary conditions
void SetBndValues1D( dTensorBC3& q)
{
    int i,m,ell;
    double tmp;
    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    
    for (ell=1; ell<=kmax; ell++)
    { 
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {        
	        // q values
	        for (m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i+melems,m,ell);
		        q.set(i,m,ell, tmp );
            }
                
        }
        // ***********************************************  

        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(melems+1); i<=(melems+mbc); i++)
        {        
	        // q values
	        for (m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i-melems,m,ell);    
		        q.set(i,m,ell, tmp );
            }
                
        }
        // ***********************************************
    }
}

/*
 *
//////////////////////////////////////////////////////////////////////////////
// Compute q_x using finite differences of moments. 
// Output is written to qout.
//
//    qout(:,1,:) = q_x
//    qout(:,2,:) = q_xx
//
//////////////////////////////////////////////////////////////////////////////
void FiniteDiff(int mestart, int meend, const dTensorBC3& q, dTensorBC3& qout)
{
    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();

    const double dx     = dogParamsCart2.get_dx();
    const double two_dx = 2.0*dx;
    const double dx2    = dx*dx;

    for (int i=1; i<=mx; i++)
    for( int me=mestart; me <= meend; me=me+2 )
    {
        dTensor1 dA ( 5  );
        dTensor1 dA2( 5  );
        dA.setall   ( 0. );
        dA2.setall  ( 0. );


        int ip1 = i+1;
        int im1 = i-1;
        int ic  = i;

        if (i<2)
        {
            im1 = im1+mx;
        }
        if (i>(mx-1))
        {
            ip1 = ip1-mx;
        }

        for (int k=1; k<=kmax; k++)
        {
            dA.set(  k, (q.get(ip1,1,k)-q.get(im1,1,k))/two_dx );
            dA2.set( k, (q.get(ip1,1,k)-2.0*q.get(ic,1,k)+q.get(im1,1,k))/dx2 );
        }

        switch(kmax)
        {
            case 5:
                qout.set(i,me,  5,  dA.get(5) );
                qout.set(i,me+1,5, dA2.get(5) );

            case 4:
                qout.set(i,me,  4,  dA.get(4) );
                qout.set(i,me+1,4, dA2.get(4) );

            case 3:
                qout.set(i,me,  3,  dA.get(3) - 14.0*sq5* dA.get(5) );
                qout.set(i,me+1,3, dA2.get(3) - 7.0*sq5*dA2.get(5) );

            case 2:
                qout.set(i,me,  2,  dA.get(2) - 10.0*sq7/sq3* dA.get(4) );
                qout.set(i,me+1,2, dA2.get(2) - 5.0*sq7/sq3*dA2.get(4) );

            case 1:
                qout.set(i,me,  1,  dA.get(1) - 2.0*sq5* dA.get(3) + 78.0* dA.get(5) );
                qout.set(i,me+1,1, dA2.get(1) - sq5*dA2.get(3) + 11.0*dA2.get(5) );

                break;
        }

    }

}

*/
