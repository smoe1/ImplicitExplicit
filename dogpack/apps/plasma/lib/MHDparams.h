#ifndef MHDPARAMS_H
#define MHDPARAMS_H
using namespace std;

class IniDocument;
// these parameters are related to PlasmaParams
struct MHDparams{
 private:
   // slowing_time^{-1} = slowing_rate
   // = resistivity * number_density * e^2 / reduced_mass
   double gamma;
   double resistivity;
   double viscosity;
   double cc; // numerical parameter
 //
 // methods
 //
 private:
   void checkParameters();
   void reportParameters();
 public:

   void init(IniDocument& ini_doc);

   double get_gamma       () const {return gamma ;}
   double get_resistivity () const {return resistivity ;}
   double get_viscosity   () const {return viscosity  ;}
   double get_cc          () const {return cc         ;}
};

extern MHDparams mhdParams;

#endif
