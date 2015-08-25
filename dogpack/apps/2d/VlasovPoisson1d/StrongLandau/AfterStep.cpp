#include<cmath>
#include "dogdefs.h"
#include "DogSolverCart2.h"

// Function that is called after each time step
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver)
{
    int i,j;
    double x,y; 
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);

    double qp, tmp, Se;
    int k;

    for( i=1-mbc; i <= mx + mbc; i++ )
    for( j=1-mbc; j <= my + mbc; j++ )
    {
        qp = pow( q.get(i,j,1,kmax), 2 );
        tmp = 0.0;
        for( k = 1; k <= kmax; k++ )
        {
            tmp += pow(q.get(i,j,1,k), 2 );
        }

        Se = qp / tmp;  // what about tmp = 0?

        aux.set(i,j,3,1,Se);
    }

}
