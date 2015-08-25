#include <cmath>
#include "dogdefs.h"
#include "AcousticParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
        const dTensor1& Q_ave, const dTensor2& Qvals, dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);
    const int kmax = Qvals.getsize(2)+1;

    // Direction
    int mu,mv;
    if (ixy==1)
    {  
        mu = 2;
        mv = 3;
    }
    else
    {
        mu = 3;
        mv = 2;
    }

    // Project onto left eigenvectors
    const double c = acousticParams.get_c();
    for (int k=1; k<=(kmax-1); k++)
    {

        Wvals.set(1,k, ( (-1.)*Qvals.get(1,k) + c * Qvals.get(mu,k) ) / (2.0*c) );
        Wvals.set(2,k, Qvals.get(mv,k) );
        Wvals.set(3,k, ( ( 1.)*Qvals.get(1,k) + c * Qvals.get(mu,k) ) / (2.0*c) );

    }

}
