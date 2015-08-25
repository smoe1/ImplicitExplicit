#ifndef TWOFLUIDPARAMS_H
#define TWOFLUIDPARAMS_H

class IniDocument;
struct TwoFluidParams
{
public:
  double gamma;
  double mass_ratio;
  double light_speed;
  void init(IniDocument& ini_doc);
};

extern TwoFluidParams twofluidParams;

#endif
