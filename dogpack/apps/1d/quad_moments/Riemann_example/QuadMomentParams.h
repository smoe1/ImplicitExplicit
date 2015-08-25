#ifndef _QUADMOMENTPARAMS_H_
#define _QUADMOMENTPARAMS_H_

class IniDocument;
class QuadMomentParams
{   
 public:
  int use_collision;
  double epsilon;


  void init(IniDocument& ini_doc);
  void write_quadMomenthelp(const char* filename);
};
extern QuadMomentParams quadMomentParams;

#endif
