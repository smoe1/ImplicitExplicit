#ifndef _INITIALPARAMS_H_
#define _INITIALPARAMS_H_

class IniDocument;
class InitialParams
{  
 public:
  // parameters for initial conditions
  double x0;
  double y0;
  double z0;
  double width;
  double u;
  double v;
  double w;
  void init(IniDocument& ini_doc);

 private:
  void write_plot_help();
};
extern InitialParams initialParams;

#endif

