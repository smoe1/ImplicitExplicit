///////////////////////////////////////////////////////////////////////////////
// Dummy function here, applications who wish to override this are required to
// do so in their own application directory.
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

void IntegratePsiBackward( double tnew, const dTensor2 &xpts, const dTensor1 &tpts, dTensor3 &psi)
{

    const int mpoints = xpts.getsize(1);
    const int meqn    = psi.getsize(2);

    for( int m=1; m <= mpoints; m++ )
    for( int me=1; me <= meqn; me++ )
    {
        for( int n=1; n <= tpts.getsize(); n++ )
        {
            psi.set(m, me, n, 0.0);
        }
    }

}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void IntegratePsiForward( double tn, const dTensor2 &xpts, const dTensor1 &tpts, dTensor3 &psi)
{

    const int mpoints = xpts.getsize(1);
    const int meqn    = psi.getsize(2);

    for( int m=1; m <= mpoints; m++ )
    for( int me=1; me <= meqn; me++ )
    {
        for( int n=1; n <= tpts.getsize(); n++ )
        {
            psi.set(m, me, n, 0.0);
        }
    }

}
///////////////////////////////////////////////////////////////////////////////
