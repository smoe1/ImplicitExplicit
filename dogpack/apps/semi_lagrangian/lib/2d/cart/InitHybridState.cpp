///////////////////////////////////////////////////////////////////////////////
// Dummy function here, applications who wish to override this are required to
// do so in their own application directory.
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParamsCart2.h"
#include "SLState.h"

void InitHybridState( 
    const dTensorBC4& q, const dTensorBC4& aux, SL_state& sl_state )
{

    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int maux = aux.getsize(2);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();

    // cell widths are uniform.
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // save node1d
    dTensor2* node1d = sl_state.node1d;
    for(int i=1; i<=(mx+1); i++)
    {
        node1d->set(i, 1, dogParamsCart2.get_xl(i)); //node.get(i,1,1) );
    }

}// end of function InitSLState
///////////////////////////////////////////////////////////////////////////////
