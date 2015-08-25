#ifndef _VLASOVPARAMS_H_
#define _VLASOVPARAMS_H_

class IniDocument;
class VlasovParams
{   
    public:

        double alpha  ;  // perturbation constant
        double temp   ;  // temperature of plasma
        double rho0   ;  // initial background density for maxwellian
        double k      ;  // phase shift added in (strong/weak/two-stream parameter)

        void init(IniDocument& ini_doc);


        // methods 
        const double& get_integration_constant(void);
        void   set_integration_constant(double c);

        const double& get_avg_momentum(void);
        void   set_avg_momemtum(double c);

    private:

        double integration_constant;  // background ion density
        double rho_u0 ;               // average momentum

};
extern VlasovParams vlasovParams;
#endif
