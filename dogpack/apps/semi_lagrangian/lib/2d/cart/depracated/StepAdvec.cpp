#include <cmath>
#include <iostream>
#include <pthread.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "discontCell.h"

///////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation:
//
//            q_t + u(y) q_x + v(x) q_y = 0 
//
// forward in time using an operator split method.  Time integration is
// only performed in one direction - x or y.
//
// Parameters
//
// dt - time step taken in each equation
// node - nodal points for grid.  method assumes uniform grid
// qold(mx,my,meqn,kmax,morder) - q before translation (last index is for SDC)
// qnew - q after being translated and projected onto the legendre basis
// advec_vel(mx, my, mpoints1d, 2) - advection velocity.
// direction - direction for velocity to take place
//			== 1: travel in x-direction (i.e. v = 0)
//			== 2: travel in y-direction (i.e. u = 0)
//
//  u1(my, mpoints1d) = u(y) - the 1d gauss point values for speeds
//  u2(mx, mpoints1d) = v(x) - the 1d gauss point values for speeds
//
///////////////////////////////////////////////////////////////////////////////


// Wrapper function to split StepAdvec into chunks
void StepAdvec(const double& dt, const dTensor3& node, 
        const dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
        const dTensorBC2& u1, const dTensorBC2& u2, int direction )
{
    void StepAdvec1(int j_start, int j_end, const double& dt, const dTensor3& node, 
            const dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
            const dTensorBC2& u1, const dTensorBC2& u2 );
    void StepAdvec2(int i_start, int i_end, const double& dt, const dTensor3& node, 
            const dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
            const dTensorBC2& u1, const dTensorBC2& u2 );

    void SetBndValues(const dTensor3& node, dTensorBC4& q, dTensorBC4& aux);

    //-local parameters -----------------------------------------------------//
    int i,j, k, me, io, jo, num_cell_shift;
    int mx   = qnew.getsize(1);
    int my   = qnew.getsize(2);
    int meqn = qnew.getsize(3);
    int kmax = qnew.getsize(4);
    int mbc  = qnew.getmbc();
    int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    int mpoints1d = int(sqrt(mpoints));

    switch( direction )
    {
        case 1:

            for(int n=1; n <= m_cuts.getsize(); n++)
            {
                StepAdvec1(1, mx, dt, node, qold, qnew, aux, u1, u2);
            }
                break;

        case 2:

            m_cuts.set(2, (int)(mx/2) );
            m_cuts.set(3, mx );
            for(int n=1; n <= m_cuts.getsize(); n++)
            {
                StepAdvec2(2, my, dt, node, qold, qnew, aux, u1, u2);
            }

                break;
    }

    SetBndValues(node, aux, qnew);

}

// this function can be run concurrently
void StepAdvec1(int j_start, int j_end, const double& dt, const dTensor3& node, 
        const dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
        const dTensorBC2& u1, const dTensorBC2& u2)
{

    //-local parameters -----------------------------------------------------//
    int i,j, k, me, io, jo, num_cell_shift;
    int mx   = qnew.getsize(1);
    int my   = qnew.getsize(2);
    int meqn = qnew.getsize(3);
    int kmax = qnew.getsize(4);
    int mbc  = qnew.getmbc();
    int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    int mpoints1d = int(sqrt(mpoints));

    // cell widths are uniform.  need only look at cell (1,1)
    double dx = node.get(2,1,1)-node.get(1,1,1);
    double dy = node.get(1,2,2)-node.get(1,1,2);
    double x1c = node.get(1,1,1) + dx/2.0;  // center of grid cell (1,1)
    double y1c = node.get(1,1,2) + dy/2.0;
    double scutx = 0.0;
    double scuty = 0.0;
    double xcut = 0.0;
    double ycut = 0.0;
    double speedx = 0.0;
    double speedy = 0.0;

    //legendre weight for "left/bottom" half of cell
    dTensor2 qleft(mpoints1d, kmax);		

    //legendre weights for "right/top" half of cell
    dTensor2 qright(mpoints1d, kmax);

    dTensor1 qcell(kmax);		//legendre weights for current cell
    dTensor1 scut(mpoints1d);   //location of discontinuity
    discontCell* cell = new discontCell(kmax, 1);
    int method1 = int((sqrt(1+8*kmax)-1)/2);
    //-----------------------------------------------------------------------//

    // Function definitions
    int iMod(int,int);			//better mod function - returns pos numbers
    double Max(double,double);

    //loop over each equation
    for(me=1; me<= meqn; me++)
    {

        //loop over every cell
        for(j=j_start; j<=j_end; j++)
        {
            for(int n=1; n<= mpoints1d; n++)
            {

                speedx = u1.get(j,n);
                num_cell_shift = (int)( floor(speedx*dt/dx) );

                //location of where discontinuity occurs in cell indexed by (1,1)
                xcut = (x1c-dx/2.0) + (-dx*num_cell_shift + speedx*dt);

                //location of discontinuity (in [-1,1])
                scutx = 2.0/dx*(xcut-x1c);
                if( fabs(scutx -1) < 1.0e-15 )
                {
                    scutx = -1.0;
                    num_cell_shift++;
                }

                scut.set(n,scutx);
            }

            cell->setScut(scut);
            for(i=1; i<=mx; i++)
            {
                for(int n=1; n<= mpoints1d; n++)
                {

                    speedx = u1.get(j,n);
                    num_cell_shift = (int)( floor(speedx*dt/dx) );

                    //location of where discontinuity occurs in cell indexed by (1,1)
                    xcut = (x1c-dx/2.0) + (-dx*num_cell_shift + speedx*dt);

                    //location of discontinuity (in [-1,1])
                    scutx = 2.0/dx*(xcut-x1c);
                    if( fabs(scutx -1) < 1.0e-15 )
                    {
                        scutx = -1.0;
                        num_cell_shift++;
                    }

                    scut.set(n,scutx);


                    /////////////////////////////////////////////////////
                    // Set ql and qr
                    /////////////////////////////////////////////////////
                    //io = 'old' cell that provides data for current cell i
                    io = (int)(i - num_cell_shift);

                    //periodic boundary conditions enforced here
                    io = iMod((io-1),mx) + 1;	

                    for(k= 1; k<=kmax; k++)
                    {

                        qright.set(n, k, qold.get(io, j, me, k) );
                        if( io - 1 == 0)
                        {
                            qleft.set(n, k, qold.get(mx, j, me, k) );
                        }else
                        {
                            qleft.set(n, k, qold.get(io-1, j, me, k) );
                        }
                    }//end of loop over each polynomial

                }

                //create a discontinuous cell and project onto legendre basis
                cell->project(qleft, qright, qcell);

                //save legendre weights into qnew
                for(k=1; k<=kmax; k++)
                {
                    qnew.set(i,j, me, k, qcell.get(k));
                }

            }
        }

    }//end of loop over each eqn	

}// end of function StepAdvecConstCoeff.cpp

void StepAdvec2(int i_start, int i_end, const double& dt, const dTensor3& node, 
        const dTensorBC4& qold, dTensorBC4& qnew, dTensorBC4& aux, 
        const dTensorBC2& u1, const dTensorBC2& u2)
{

    //-local parameters -----------------------------------------------------//
    int i,j, k, me, io, jo, num_cell_shift;
    int mx   = qnew.getsize(1);
    int my   = qnew.getsize(2);
    int meqn = qnew.getsize(3);
    int kmax = qnew.getsize(4);
    int mbc  = qnew.getmbc();
    int mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    int mpoints1d = int(sqrt(mpoints));

    // cell widths are uniform.  need only look at cell (1,1)
    double dx = node.get(2,1,1)-node.get(1,1,1);
    double dy = node.get(1,2,2)-node.get(1,1,2);
    double x1c = node.get(1,1,1) + dx/2.0;  // center of grid cell (1,1)
    double y1c = node.get(1,1,2) + dy/2.0;
    double scutx = 0.0;
    double scuty = 0.0;
    double xcut = 0.0;
    double ycut = 0.0;
    double speedx = 0.0;
    double speedy = 0.0;

    //legendre weight for "left/bottom" half of cell
    dTensor2 qleft(mpoints1d, kmax);		

    //legendre weights for "right/top" half of cell
    dTensor2 qright(mpoints1d, kmax);

    dTensor1 qcell(kmax);		//legendre weights for current cell
    dTensor1 scut(mpoints1d);   //location of discontinuity
    discontCell* cell = new discontCell(kmax, 2);
    int method1 = int((sqrt(1+8*kmax)-1)/2);
    //-----------------------------------------------------------------------//

    // Function definitions
    int iMod(int,int);			//better mod function - returns pos numbers
    double Max(double,double);

    //loop over each equation
    for(me=1; me<= meqn; me++)
    {

        //loop over every cell
        for(i=i_start; i<=i_end; i++)
        {
            for( int n=1; n<= mpoints1d; n++)
            {
                speedy = u2.get(i, n);
                num_cell_shift = (int) ( floor(speedy * dt/ dy) );

                //location of where discontinuity occurs in cell indexed by (1,1)
                ycut = (y1c-dy/2.0) + (-dy*num_cell_shift + speedy*dt);

                //location of discontinuity (in [-1,1])
                scuty = 2.0/dy*(ycut-y1c);

                if( fabs(scuty -1) < 1.0e-14 )
                {
                    scuty = -1.0;
                    num_cell_shift++;
                }

                scut.set(n,scuty);

            }

            // this is the same for every collumn, and is an expensive
            // call
            cell->setScut(scut);

            for(j=1; j<=my; j++)
            {

                for( int n=1; n<= mpoints1d; n++)
                {

                    /////////////////////////////////////////////////////
                    // Set ql and qr
                    /////////////////////////////////////////////////////					
                    //jo = 'old' cell that provides data for current cell j
                    speedy = u2.get(i, n);
                    num_cell_shift = (int) ( floor(speedy * dt/ dy) );
                    jo = (int)(j - num_cell_shift);

                    //save legendre weight from qold using
                    //using periodic boundary conditions
                    jo = iMod((jo-1),my) + 1;

                    for(k= 1; k<=kmax; k++)
                    {

                        qright.set(n, k, qold.get(i, jo, me, k) );
                        if( jo - 1 == 0)
                        {
                            qleft.set(n, k, qold.get(i, my, me, k) );
                        }else
                        {
                            qleft.set(n, k, qold.get(i, jo-1, me, k) );
                        }
                    }//end of loop over each polynomial
                }


                //create a discontinuous cell and project onto legendre basis
                cell->project(qleft, qright, qcell);

                //save legendre weights into qnew
                for(k=1; k<=kmax; k++)
                {
                    qnew.set(i,j, me, k, qcell.get(k));
                }
            }
        }//end of loop over each cell ( advection in y-direction)

        delete cell;

    }//end of loop over each eqn	

}// end of function StepAdvecConstCoeff.cpp
