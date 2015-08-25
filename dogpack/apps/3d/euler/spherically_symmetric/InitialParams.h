#ifndef _INITIALPARAMS_H_
#define _INITIALPARAMS_H_

class IniDocument;
class InitialParams
{  
 public:
  // parameters for initial conditions
  double rhol, u1l, u2l, u3l, pl;
  double rhor, u1r, u2r, u3r, pr;

  // derived parameters
  double m1l, m2l, m3l, El;
  double m1r, m2r, m3r, Er;

  // functions
  void init(IniDocument& ini_doc);
  void set_derived_params(double gamma);

 private:
  void write_plot_help();
};
extern InitialParams initialParams;

#endif

