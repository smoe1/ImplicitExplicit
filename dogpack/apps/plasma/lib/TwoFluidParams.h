#ifndef TwoFluidParams_h
#define TwoFluidParams_h

struct TwoFluidParams{
    double mass_ratio;
    double ion_mass;
    double elc_mass;
  public:
    TwoFluidParams()
    {
        mass_ratio=1.;
        ion_mass=1.;
        elc_mass=1.;
    }
};

#endif
