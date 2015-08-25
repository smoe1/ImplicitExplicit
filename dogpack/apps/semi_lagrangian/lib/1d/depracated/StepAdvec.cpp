#include <cmath>
#include "dog_math.h"
#include "tensors.h"

//////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation forward in time.
//
// Parameters:
//
// dt - time step taken in each equation
// node - list of node points: node(1,1) = x_low, ... node(mx+1,1) = x_high
// auxvals - advection speed in each equation
// smax - maximum wave speed in each cell -- this is constant for advection
// eqution!
// qvals - q before translation
// qnew - the returned values (q(i,m,k) = Q^{(k)}_i for the mth equation.
//        this is q after being translated and projected onto the legendre
//        after we use two halves to compute the integration
//////////////////////////////////////////////////////////////////////////////
void StepAdvec(const double& dt, const dTensor2& node,
        dTensorBC3& auxvals, dTensorBC1& smax,
        const dTensorBC3& qvals, dTensorBC3& qnew)
{
    //-function declarations-------------------------------------------------//
    void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
    //-----------------------------------------------------------------------//

    //-local parameters -----------------------------------------------------//
    const int melems  = qvals.getsize(1);     //number of elements in grid
    const int meqn    = qvals.getsize(2);     //number of equations
    const int morder = qvals.getsize(3);      //order of accuracy in space
    const int mbc = qvals.getmbc();        //number of ghost cells
    //-----------------------------------------------------------------------//

    const int kmax = qnew.getsize(3);

    const double dx = node.get(2,1)-node.get(1,1);
    const double speed = auxvals.get(1,1,1);

    //loop over cell and equation
#pragma omp parallel for
    for ( int i=1; i<=melems; i++)
        for ( int me=1; me<= meqn; me++)
        {
            {

                dTensor1 qleft(morder);        //legendre weights on grid cell j-1
                dTensor1 qright(morder);       //legendre weights on grid cell j

                // save max wave speed == fabs(speed)  -- this is a constant!
                smax.set(i, Max(smax.get(i), fabs(speed) ) );

                double Ql1, Ql2, Ql3, Ql4, Ql5;
                double Q1, Q2, Q3, Q4, Q5;

                //find cell that contributes to cell xi
                int numCellShift = (int)( floor(speed*dt/dx) );

                // periodicity is enforced here!
                int j = (int)(i - numCellShift);
                j = iMod((j-1),melems) + 1;

                //save legendre weights from qold
                for(int k = 1; k<=morder; k++)
                {
                    qright.set( k, qvals.get(j, me, k) );
                    if(j-1 == 0)
                    {
                        qleft.set(k, qvals.get(melems,me,k) );
                    }
                    else
                    {
                        qleft.set(k, qvals.get(j-1,me,k) );
                    }
                }//end of loop over each polynomial

                double nu = speed * dt / dx - numCellShift;
                double tmp = 0.0;
                switch( kmax )
                {

                    case 5:
                        Ql2 = qleft. get(2);   
                        Ql1 = qleft. get(1);
                        Q2  = qright.get(2);   
                        Q1  = qright.get(1);
                        Ql3 = qleft. get(3);   
                        Q3  = qright.get(3);
                        Ql4 = qleft. get(4);   
                        Q4  = qright.get(4);
                        Ql5 = qleft. get(5);   
                        Q5  = qright.get(5);

                        // k == 1
                        tmp = (double) 
                            (42*Ql5-42*Q5)*pow(nu,5)+(5*Q4*sqrt(7)-5*Ql4*sqrt(7)+105*Q5-105*Ql5)*pow(nu,4)+(-90*Q5+10*Ql4*sqrt(7)+90*Ql5-10*Q4*sqrt(7)+2*Ql3*sqrt(5)-2*Q3*sqrt(5))*pow(nu,3)+(-3*Ql3*sqrt(5)-6*Ql4*sqrt(7)+Q2*sqrt(3)+30*Q5-Ql2*sqrt(3)-30*Ql5+6*Q4*sqrt(7)+3*Q3*sqrt(5))*pow(nu,2)+(3*Ql5+Ql3*sqrt(5)+Ql4*sqrt(7)+Ql2*sqrt(3)-Q1-Q2*sqrt(3)+Ql1-Q4*sqrt(7)-3*Q5-Q3*sqrt(5))*nu+Q1;
                        qnew.set(i,me,1, tmp);

                        // k == 2
                        tmp = (double)
                            (-14*sqrt(3)*Q5+14*sqrt(3)*Ql5)*pow(nu,6)+(2*sqrt(3)*Q4*sqrt(7)-84*sqrt(3)*Ql5-2*sqrt(3)*Ql4*sqrt(7))*pow(nu,5)+(-sqrt(3)*Q3*sqrt(5)+sqrt(3)*Ql3*sqrt(5)+60*sqrt(3)*Q5+10*sqrt(3)*Ql4*sqrt(7)+150*sqrt(3)*Ql5)*pow(nu,4)+(-4*sqrt(3)*Ql3*sqrt(5)-110*sqrt(3)*Ql5+2*Q2-14*sqrt(3)*Ql4*sqrt(7)-6*sqrt(3)*Q4*sqrt(7)-2*Ql2-70*sqrt(3)*Q5)*pow(nu,3)+(6*Ql2-sqrt(3)*Q1+4*sqrt(3)*Ql3*sqrt(5)+33*sqrt(3)*Ql5+7*sqrt(3)*Ql4*sqrt(7)+27*sqrt(3)*Q5+sqrt(3)*Ql1+5*sqrt(3)*Q4*sqrt(7)+2*sqrt(3)*Q3*sqrt(5))*pow(nu,2)+(sqrt(3)*Q1-3*sqrt(3)*Q5-sqrt(3)*Ql3*sqrt(5)-sqrt(3)*Ql1-sqrt(3)*Q3*sqrt(5)-sqrt(3)*Ql4*sqrt(7)-3*sqrt(3)*Ql5-3*Ql2-3*Q2-sqrt(3)*Q4*sqrt(7))*nu+Q2;
                        qnew.set(i,me,2, tmp );

                        // k == 3 
                        tmp = (double)
                            (-12*sqrt(5)*Q5+12*sqrt(5)*Ql5)*pow(nu,7)+(-84*sqrt(5)*Ql5+2*Q4*sqrt(7)*sqrt(5)-2*Ql4*sqrt(7)*sqrt(5))*pow(nu,6)+(-6*Q3+30*sqrt(5)*Q5+222*sqrt(5)*Ql5+6*Ql3+12*Ql4*sqrt(7)*sqrt(5))*pow(nu,5)+(-26*Ql4*sqrt(7)*sqrt(5)-4*Q4*sqrt(7)*sqrt(5)-30*Ql3-sqrt(5)*Ql2*sqrt(3)-270*sqrt(5)*Ql5+sqrt(5)*Q2*sqrt(3))*pow(nu,4)+(2*sqrt(5)*Ql1-2*sqrt(5)*Q1+4*sqrt(5)*Ql2*sqrt(3)-36*sqrt(5)*Q5+24*Ql4*sqrt(7)*sqrt(5)+156*sqrt(5)*Ql5+10*Q3+50*Ql3)*pow(nu,3)+(3*sqrt(5)*Q1-3*sqrt(5)*Ql1-39*sqrt(5)*Ql5-4*sqrt(5)*Ql2*sqrt(3)+21*sqrt(5)*Q5+3*Q4*sqrt(7)*sqrt(5)-30*Ql3-2*sqrt(5)*Q2*sqrt(3)-9*Ql4*sqrt(7)*sqrt(5))*pow(nu,2)+(sqrt(5)*Ql1+sqrt(5)*Ql2*sqrt(3)+sqrt(5)*Q2*sqrt(3)-sqrt(5)*Q1+Ql4*sqrt(7)*sqrt(5)-5*Q3+5*Ql3-3*sqrt(5)*Q5+3*sqrt(5)*Ql5-Q4*sqrt(7)*sqrt(5))*nu+Q3;
                        qnew.set(i,me,3, tmp);

                        // k == 4 
                        tmp = (double)
                            (-15*sqrt(7)*Q5+15*sqrt(7)*Ql5)*pow(nu,8)+(-120*sqrt(7)*Ql5-20*Ql4+20*Q4)*pow(nu,7)+(2*sqrt(7)*Ql3*sqrt(5)-2*sqrt(7)*Q3*sqrt(5)+140*Ql4+36*sqrt(7)*Q5+384*sqrt(7)*Ql5)*pow(nu,6)+(-624*sqrt(7)*Ql5-12*sqrt(7)*Ql3*sqrt(5)-42*Q4-378*Ql4-2*sqrt(7)*Ql2*sqrt(3)+2*sqrt(7)*Q2*sqrt(3))*pow(nu,5)+(10*sqrt(7)*Ql2*sqrt(3)+5*sqrt(7)*Ql1-5*sqrt(7)*Q1+4*sqrt(7)*Q3*sqrt(5)+490*Ql4+26*sqrt(7)*Ql3*sqrt(5)-30*sqrt(7)*Q5+540*sqrt(7)*Ql5)*pow(nu,4)+(-308*Ql4-240*sqrt(7)*Ql5-6*sqrt(7)*Q2*sqrt(3)-14*sqrt(7)*Ql2*sqrt(3)-10*sqrt(7)*Ql1-24*sqrt(7)*Ql3*sqrt(5)+28*Q4+10*sqrt(7)*Q1)*pow(nu,3)+(9*sqrt(7)*Ql3*sqrt(5)+84*Ql4+48*sqrt(7)*Ql5+12*sqrt(7)*Q5+7*sqrt(7)*Ql2*sqrt(3)+5*sqrt(7)*Q2*sqrt(3)-3*sqrt(7)*Q3*sqrt(5)+6*sqrt(7)*Ql1-6*sqrt(7)*Q1)*pow(nu,2)+(-sqrt(7)*Ql1-sqrt(7)*Ql2*sqrt(3)-3*sqrt(7)*Ql5+sqrt(7)*Q3*sqrt(5)-7*Q4-sqrt(7)*Q2*sqrt(3)-7*Ql4-sqrt(7)*Ql3*sqrt(5)-3*sqrt(7)*Q5+sqrt(7)*Q1)*nu+Q4;
                        qnew.set(i,me,4, tmp);

                        // k == 5 
                        tmp = (double)
                            (-70*Q5+70*Ql5)*pow(nu,9)+(-15*Ql4*sqrt(7)-630*Ql5+15*Q4*sqrt(7))*pow(nu,8)+(-12*Q3*sqrt(5)+180*Q5+120*Ql4*sqrt(7)+12*Ql3*sqrt(5)+2340*Ql5)*pow(nu,7)+(-4620*Ql5-384*Ql4*sqrt(7)-14*Ql2*sqrt(3)-84*Ql3*sqrt(5)-36*Q4*sqrt(7)+14*Q2*sqrt(3))*pow(nu,6)+(624*Ql4*sqrt(7)+5202*Ql5+42*Ql1-162*Q5+222*Ql3*sqrt(5)+30*Q3*sqrt(5)-42*Q1+84*Ql2*sqrt(3))*pow(nu,5)+(-150*Ql2*sqrt(3)-270*Ql3*sqrt(5)-3330*Ql5+105*Q1+30*Q4*sqrt(7)-60*Q2*sqrt(3)-540*Ql4*sqrt(7)-105*Ql1)*pow(nu,4)+(240*Ql4*sqrt(7)+1140*Ql5+70*Q2*sqrt(3)-36*Q3*sqrt(5)+156*Ql3*sqrt(5)+60*Q5-90*Q1+110*Ql2*sqrt(3)+90*Ql1)*pow(nu,3)+(-12*Q4*sqrt(7)+21*Q3*sqrt(5)-30*Ql1-33*Ql2*sqrt(3)-48*Ql4*sqrt(7)+30*Q1-27*Q2*sqrt(3)-180*Ql5-39*Ql3*sqrt(5))*pow(nu,2)+(9*Ql5+3*Ql3*sqrt(5)+3*Ql4*sqrt(7)+3*Ql2*sqrt(3)+3*Q2*sqrt(3)+3*Ql1-9*Q5+3*Q4*sqrt(7)-3*Q1-3*Q3*sqrt(5))*nu+Q5;
                        qnew.set(i,me,5, tmp);


                        break;

                    case 4:
                        Ql2 = qleft. get(2);   
                        Ql1 = qleft. get(1);
                        Q2  = qright.get(2);   
                        Q1  = qright.get(1);
                        Ql3 = qleft. get(3);   
                        Q3  = qright.get(3);
                        Ql4 = qleft. get(4);   
                        Q4  = qright.get(4);

                        // k == 1
                        tmp = (double) 
                            (5*Q4*sqrt(7)-5*Ql4*sqrt(7))*pow(nu,4)+(-10*Q4*sqrt(7)-2*Q3*sqrt(5)+10*Ql4*sqrt(7)+2*Ql3*sqrt(5))*pow(nu,3)+(6*Q4*sqrt(7)+3*Q3*sqrt(5)-Ql2*sqrt(3)-6*Ql4*sqrt(7)-3*Ql3*sqrt(5)+Q2*sqrt(3))*pow(nu,2)+(Ql4*sqrt(7)+Ql1-Q2*sqrt(3)+Ql3*sqrt(5)-Q4*sqrt(7)+Ql2*sqrt(3)-Q1-Q3*sqrt(5))*nu+Q1;
                        qnew.set(i,me,1, tmp);

                        // k == 2
                        tmp = (double)
                            (-2*sqrt(3)*Ql4*sqrt(7)+2*sqrt(3)*Q4*sqrt(7))*pow(nu,5)+(-sqrt(3)*Q3*sqrt(5)+10*sqrt(3)*Ql4*sqrt(7)+sqrt(3)*Ql3*sqrt(5))*pow(nu,4)+(-14*sqrt(3)*Ql4*sqrt(7)+2*Q2-6*sqrt(3)*Q4*sqrt(7)-2*Ql2-4*sqrt(3)*Ql3*sqrt(5))*pow(nu,3)+(sqrt(3)*Ql1+7*sqrt(3)*Ql4*sqrt(7)+6*Ql2-sqrt(3)*Q1+5*sqrt(3)*Q4*sqrt(7)+2*sqrt(3)*Q3*sqrt(5)+4*sqrt(3)*Ql3*sqrt(5))*pow(nu,2)+(-3*Q2-3*Ql2-sqrt(3)*Ql3*sqrt(5)-sqrt(3)*Q4*sqrt(7)-sqrt(3)*Ql4*sqrt(7)-sqrt(3)*Q3*sqrt(5)-sqrt(3)*Ql1+sqrt(3)*Q1)*nu+Q2;
                        qnew.set(i,me,2, tmp );

                        // k == 3 
                        tmp = (double)
                            (-2*sqrt(5)*Ql4*sqrt(7)+2*sqrt(5)*Q4*sqrt(7))*pow(nu,6)+(6*Ql3-6*Q3+12*sqrt(5)*Ql4*sqrt(7))*pow(nu,5)+(-4*sqrt(5)*Q4*sqrt(7)-sqrt(5)*Ql2*sqrt(3)+sqrt(5)*Q2*sqrt(3)-26*sqrt(5)*Ql4*sqrt(7)-30*Ql3)*pow(nu,4)+(50*Ql3+10*Q3+4*sqrt(5)*Ql2*sqrt(3)+24*sqrt(5)*Ql4*sqrt(7)-2*sqrt(5)*Q1+2*sqrt(5)*Ql1)*pow(nu,3)+(-4*sqrt(5)*Ql2*sqrt(3)-9*sqrt(5)*Ql4*sqrt(7)-2*sqrt(5)*Q2*sqrt(3)+3*sqrt(5)*Q1-3*sqrt(5)*Ql1+3*sqrt(5)*Q4*sqrt(7)-30*Ql3)*pow(nu,2)+(sqrt(5)*Ql2*sqrt(3)+sqrt(5)*Q2*sqrt(3)+sqrt(5)*Ql1+5*Ql3+sqrt(5)*Ql4*sqrt(7)-sqrt(5)*Q1-5*Q3-sqrt(5)*Q4*sqrt(7))*nu+Q3;
                        qnew.set(i,me,3, tmp);
                        // k == 4 
                        tmp = (double)
                            (-20*Ql4+20*Q4)*pow(nu,7)+(140*Ql4+2*sqrt(7)*Ql3*sqrt(5)-2*sqrt(7)*Q3*sqrt(5))*pow(nu,6)+(-378*Ql4-12*sqrt(7)*Ql3*sqrt(5)+2*sqrt(7)*Q2*sqrt(3)-2*sqrt(7)*Ql2*sqrt(3)-42*Q4)*pow(nu,5)+(4*sqrt(7)*Q3*sqrt(5)+26*sqrt(7)*Ql3*sqrt(5)+10*sqrt(7)*Ql2*sqrt(3)-5*sqrt(7)*Q1+5*sqrt(7)*Ql1+490*Ql4)*pow(nu,4)+(-6*sqrt(7)*Q2*sqrt(3)-10*sqrt(7)*Ql1+28*Q4-308*Ql4-24*sqrt(7)*Ql3*sqrt(5)-14*sqrt(7)*Ql2*sqrt(3)+10*sqrt(7)*Q1)*pow(nu,3)+(7*sqrt(7)*Ql2*sqrt(3)-3*sqrt(7)*Q3*sqrt(5)+84*Ql4+5*sqrt(7)*Q2*sqrt(3)+6*sqrt(7)*Ql1-6*sqrt(7)*Q1+9*sqrt(7)*Ql3*sqrt(5))*pow(nu,2)+(-7*Ql4+sqrt(7)*Q1-sqrt(7)*Ql2*sqrt(3)-7*Q4-sqrt(7)*Ql3*sqrt(5)-sqrt(7)*Ql1+sqrt(7)*Q3*sqrt(5)-sqrt(7)*Q2*sqrt(3))*nu+Q4;
                        qnew.set(i,me,4, tmp);

                        break;

                    case 3:
                        Ql2 = qleft. get(2);   
                        Ql1 = qleft. get(1);
                        Q2  = qright.get(2);   
                        Q1  = qright.get(1);
                        Ql3 = qleft. get(3);   
                        Q3  = qright.get(3);

                        // k == 1
                        tmp = (double) (2*Ql3*sqrt(5)-2*Q3*sqrt(5))*pow(nu,3)+(-3*Ql3*sqrt(5)+3*Q3*sqrt(5)+sqrt(3)*Q2-Ql2*sqrt(3))*pow(nu,2)+(-Q3*sqrt(5)+Ql1+Ql2*sqrt(3)+Ql3*sqrt(5)-Q1-sqrt(3)*Q2)*nu+Q1;
                        qnew.set(i,me,1, tmp);

                        // k == 2
                        tmp = (double)
                            (sqrt(3)*Ql3*sqrt(5)-sqrt(3)*Q3*sqrt(5))*pow(nu,4)+(-2*Ql2+2*Q2-4*sqrt(3)*Ql3*sqrt(5))*pow(nu,3)+(sqrt(3)*Ql1-sqrt(3)*Q1+4*sqrt(3)*Ql3*sqrt(5)+6*Ql2+2*sqrt(3)*Q3*sqrt(5))*pow(nu,2)+(-3*Ql2-sqrt(3)*Q3*sqrt(5)-3*Q2-sqrt(3)*Ql1+sqrt(3)*Q1-sqrt(3)*Ql3*sqrt(5))*nu+Q2;
                        qnew.set(i,me,2, tmp );

                        // k == 3 
                        tmp = (double)
                            (-6*Q3+6*Ql3)*pow(nu,5)+(sqrt(3)*Q2*sqrt(5)-sqrt(5)*Ql2*sqrt(3)-30*Ql3)*pow(nu,4)+(4*sqrt(5)*Ql2*sqrt(3)-2*sqrt(5)*Q1+10*Q3+2*sqrt(5)*Ql1+50*Ql3)*pow(nu,3)+(-2*sqrt(3)*Q2*sqrt(5)-4*sqrt(5)*Ql2*sqrt(3)-3*sqrt(5)*Ql1-30*Ql3+3*sqrt(5)*Q1)*pow(nu,2)+(-5*Q3-sqrt(5)*Q1+sqrt(5)*Ql2*sqrt(3)+sqrt(5)*Ql1+5*Ql3+sqrt(3)*Q2*sqrt(5))*nu+Q3;
                        qnew.set(i,me,3, tmp);
                        break;

                    case 2:

                        Ql2 = qleft.get(2);   
                        Ql1 = qleft.get(1);
                        Q2  = qright.get(2);   
                        Q1  = qright.get(1);

                        tmp = (-sqrt(3.0)*nu+sqrt(3.0)*pow(nu,2))*Q2
                            +Ql2*sqrt(3.0)*nu-Ql2*sqrt(3.0)*pow(nu,2)+nu*Ql1+Q1-Q1*nu;
                        qnew.set(i,me,1, tmp );

                        tmp = (sqrt(3.0)*nu-sqrt(3.0)*pow(nu,2))*Q1-3.0*Ql2*nu+6.0*Ql2*pow(nu,2)
                            -3.0*Q2*nu-sqrt(3.0)*Ql1*nu+sqrt(3.0)*Ql1*pow(nu,2)
                            +Q2-2*Ql2*pow(nu,3)+2*Q2*pow(nu,3);
                        qnew.set(i,me,2,tmp);
                        break;

                    case 1:
                        qnew.set(i,me,1, nu * qleft.get(1) + (1.0-nu) * qright.get(1) );

                        break;

                }

            }//end of loop over each cell
        }//end of loop over each eqn  

    SetBndValues(node, auxvals, qnew);

}
