#ifndef _INITIALPARAMS_H_
#define _INITIALPARAMS_H_

class IniDocument;
class InitialParams
{  
    public:
        // parameters for initial conditions
        double x0,y0,z0,w0;
        double width;
        double u1,u2,u3,u4;
        void init(IniDocument& ini_doc);

    private:
        void write_plot_help();
};
extern InitialParams initialParams;

#endif

