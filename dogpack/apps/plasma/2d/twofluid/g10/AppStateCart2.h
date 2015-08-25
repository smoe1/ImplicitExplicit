#ifndef AppStateCart2_h
#define AppStateCart2_h
#include "DogStateCart2.h"

class AppSolverCart2;
class AppStateCart2 : public DogStateCart2
{
 public:
  virtual const char* get_classname(){return (const char*) "AppStateCart2";}
  AppStateCart2(){};
 private:
  AppStateCart2& operator=(const AppStateCart2& in); // disabled
 private:
  AppStateCart2(const AppStateCart2&in, CopyMode::Enum copyMode):
    DogStateCart2(in, copyMode) { }
  virtual DogState* clone(CopyMode::Enum copyMode)const
  { return new AppStateCart2(*this, copyMode); }
  //virtual void AfterFullTimeStep(double dt);
  virtual int advanceSplitTimeStep(double dt);
  virtual void BeforeStage(double dt);
  virtual bool AfterAdvanceStage();
  virtual void ReportAfterStep()const;

  // To grant access to the solver we can simply put a pointer
  // to the solver in the DogState class.  This avoids forcing
  // the solver to be singleton.  Again, we cast to avoid
  // the need to propagate pointers.
  //
  const AppSolverCart2& get_solver()const
  { return (const AppSolverCart2&) DogState::get_solver();}
  AppSolverCart2& fetch_solver()
  { return (AppSolverCart2&) DogState::fetch_solver();}
};

#endif // AppStateCart2_h
