#ifndef _HYBRIDPARAMS_H_
#define _HYBRIDPARAMS_H_

class IniDocument;
class SLParams 
{   
    public:

        void init(IniDocument& ini_doc);

        // TODO - would be nice to make these strings ...
        //char* time_splitting_method;  // flag to indicate if using hybrid method

        // flag to indicate if using rk-nystrom (default =0, no )
        int rk_nystrom_flag;  

        // flag to indicate if using hybrid method - default = no
        int hybrid_flag;      
        double Rk_Cfl ;       // Runge-Kutta CFL number - used for hybrid scheme 

        // methods 
//      const double& get_integration_constant(void);
//      void   set_integration_constant(double c);

//  private:

};
extern SLParams slParams;
#endif
