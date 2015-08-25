
struct AppParams{

   // directly configured
   //
   double gamma;
   double mass_ratio;
   double debye;
   double cs_light;
   double larmor_radius;
   double ion_rlx_period;
   double elc_rlx_period;
   //
   // derived
   //
   double ion_mass;
   double elc_mass;
};

extern AppParams appParams;

