#ifndef _COLLISIONPARAMS_H_
#define _COLLISIONPARAMS_H_

class IniDocument;
class CollisionParams
{   
 public:
  double epsilon;
  void init(IniDocument& ini_doc);
  void write_collisionhelp(const char* filename);
};
extern CollisionParams collisionParams;

#endif
