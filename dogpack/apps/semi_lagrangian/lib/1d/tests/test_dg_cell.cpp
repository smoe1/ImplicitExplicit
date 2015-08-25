#include <stdio.h>

#ifndef CHECK_BOUNDS
#define CHECK_BOUNDS
#endif

#ifndef NDIMS
#define NDIMS 0
#endif

#include <cmath>
#include "tensors.h"
#include "assert.h"

#include "../dg_cell.h"

int err_code( int k )
{
    printf("  failed test %d\n", k );
    return -1;
}

int main()
{

    const double eps = 1e-14;

    dg_cell* dgc = new dg_cell();
    dTensor1 phik( 5 );

    // TEST 1:
    // check that integral over whole domain is \delta_{1, k} //
    dgc->integratePhi( -1.0, 1.0,  phik );
    assert( fabs( phik.get(1) / 2.0 - 1.0 ) < eps );
    for( int k=2; k <= 5; k++ )
    { assert( fabs( phik.get(k) ) < eps ); }
    printf("   Passed test 1\n");

    // TEST 2:
    // check interval [-1,0] //
    dgc->integratePhi( -1.0, 0.0, phik );
    assert( fabs( phik.get(1) - 1.0 ) < eps );
    assert( fabs( phik.get(2) + sqrt(3.0) / 2.0 ) < eps );
    assert( fabs( phik.get(3) - 0.0 ) < eps );
    assert( fabs( phik.get(4) - sqrt(7.0)/8.0 ) < eps );
    assert( fabs( phik.get(5) ) < eps );
    printf("   Passed test 2\n");

    // TEST 3:
    // check interval [-0.1, 0.5] //
    dgc->integratePhi( -0.1, 0.5, phik );
    assert( fabs( phik.get(1) - 3.0/5.0 ) < eps );
    assert( fabs( phik.get(2) - 3.0 * sqrt(3.0) / 25.0 < eps ) );
    assert( fabs( phik.get(3) - ( -0.237*sqrt(5.0) ) < eps ) );
    assert( fabs( phik.get(4) - ( -0.141*sqrt(7.0) ) < eps ) );
    assert( fabs( phik.get(5) - ( 113823.0 / 10000.0 ) < eps ) );
    printf("   Passed test 3\n");

    printf("   Passed 3 of 3 integration tests! \n");

    delete dgc;
    return 0;
    
}
