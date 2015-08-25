///////////////////////////////////////////////////////////////////////////////
// Dummy function here, applications who wish to override this are required to
// do so in their own application directory.
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParamsCart2.h"
#include "SLState.h"

void SetAdvecSpeedExtra( 
    const dTensorBC4& q, const dTensorBC4& aux, dTensorBC3& aux_extra, 
    SL_state& sl_state )
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux_extra.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();

    // cell widths are uniform.
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

}// end of function SetAdvecSpeedExtra
///////////////////////////////////////////////////////////////////////////////
