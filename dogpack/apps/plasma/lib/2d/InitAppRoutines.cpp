#include "assert.h"
#include "debug.h"
#include "constants.h"

// dogpack headers
#include "DogParamsCart2.h"

// plasma/lib headers
#include "PlasmaParams.h"

// 2d/plasma/lib headers
//#include "Boundaries.h"
#include "GEMparams.h"

using namespace std;

class IniDocument;
void InitPlasma2(IniDocument& ini_doc)
{
    plasmaParams.init(ini_doc);
    gemParams.init(ini_doc);

    // check parameters
    //
    bool test_equal(double a, double b);
    switch(plasmaParams.get_model())
    {
        case PlasmaModelType::p05:
        case PlasmaModelType::p10:
        case PlasmaModelType::p20:
            if(!test_equal(gemParams.get_B_guide(),0.))
            {
              eprintf("guide field must be zero for symmetric pair plasma");
            }
    }
    if(!test_equal(gemParams.get_temp_ratio()*gemParams.get_temp_ratio(),
                plasmaParams.get_mass_ratio()))
    {
        printf(
                "warning: for the GEM problem the temperature ratio"
                " should be the square root of the mass ratio:"
                "\n  temperature ratio = %e"
                "\n  mass ratio = %e\n",
                gemParams.get_temp_ratio(),plasmaParams.get_mass_ratio());
    }
}

bool resetDogParamsCart2(double domain_scaling, 
        int enforced_symmetry)
{
    bool changed = false;

    // not necessarily true for generalized problem
    //
    const double xhigh = dogParamsCart2.get_xhigh();
    const double yhigh = dogParamsCart2.get_yhigh();
    if(domain_scaling != 1.)
    {
        assert_almost_eq(xhigh, 4.*pi*domain_scaling);
        assert_almost_eq(yhigh, 2.*pi*domain_scaling);
    }

    double xlow = -xhigh;
    double ylow = -yhigh;
    // enforce horizontal domain symmetry ( X or XY )
    if((enforced_symmetry & 1) == 1)
    {
        if(xlow != 0) changed=true;
        xlow = 0;
    }
    // enforce vertical domain symmetry ( Y or XY )
    if((enforced_symmetry & 2) == 2)
    {
        if(ylow != 0) changed=true;
        ylow = 0;
    }
    // bit 4 means mx and my describe actual computational mesh
    // enforce horizontal domain symmetry ( X or XY )
    if((enforced_symmetry & 5) == 1)
    {
        // currently only support even mx
        assert((dogParamsCart2.get_mx()&1)==0);
        dogParamsCart2.reset_mx(dogParamsCart2.get_mx()/2);
        changed=true;
    }
    // enforce vertical domain symmetry ( Y or XY )
    if((enforced_symmetry & 6) == 2)
    {
        // currently only support even my
        assert((dogParamsCart2.get_my()&1)==0);
        dogParamsCart2.reset_my(dogParamsCart2.get_my()/2);
        changed=true;
    }
    dogParamsCart2.set_xlims(xlow,xhigh);
    dogParamsCart2.set_ylims(ylow,yhigh);
    return changed;
}

