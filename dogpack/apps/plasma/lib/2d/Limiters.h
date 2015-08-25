#ifndef Limiters_h
#define Limiters_h

namespace PositivityReturnType
{
  enum Enum {
    NEG_DENSITY=-2,
    NEG_PRESSURE=-1,
    UNDEFINED=0,
    ENFORCED_DENSITY=1,
    ENFORCED_PRESSURE=2,
    ENFORCED_POSITIVITY=3,
    EVALUATIONS_CONFIRMED_POSITIVITY=4,
    INTERVALS_CONFIRMED_POSITIVITY=5,
    MAX=6,
  };
};

const char* get_string(PositivityReturnType::Enum in);

// a way to define read-only constants without polluting
// the global namespace
namespace LimiterParams
{
  // These are presumed to be a larger order of magnitude
  // than EPSILON.  The corresponding primitive variables
  // are not allowed to fall below these thresholds, which
  // can incur unphysicality (e.g. loss of conservation)
  // of the order of these epsilon values.
  // 
  inline double get_rho_epsilon() { return 1e-12; }
  inline double get_pressure_epsilon() { return 1e-12; }
  // a little bigger so that that we don't have to completely
  // squash devations
  inline double get_avg_rho_epsilon()
  { return 1.25 * get_rho_epsilon(); }
  inline double get_avg_pressure_epsilon()
  { return 1.25 * get_pressure_epsilon(); }
  // use a smaller value when verifying so we don't get a violation
  // due to numerical error
  //inline double get_verify_rho_epsilon()
  //{ return       .9 * get_rho_epsilon(); }
  //inline double get_verify_avg_rho_epsilon()
  //{ return       .9 * get_avg_rho_epsilon(); }
  inline double get_verify_pressure_epsilon()
  { return       .9 * get_pressure_epsilon(); }
  inline double get_verify_avg_pressure_epsilon()
  { return       .9 * get_avg_pressure_epsilon(); }
  
  inline double get_pos3_shift(){return 1e-19;}
}

class IntervalArray;
class PositivityPointInfo;
class CellPositivityLimiter
{
 protected: // too lazy to write accessors
  int i;
  int j;
  int moffset;
  dTensorBC4& q;
  const PositivityPointInfo& points;
 private:
  //double get_minimum_density_at_coarse_points()const;
  double get_minimum_density_at_points(const dTensor2& phi)const;
 public:
  void get_Q_at_point(dTensor1& Q, int mp)const;
  void get_dQ_at_positivity_points(dTensor2& dQ)const;
  void get_Q0_in_cell(dTensor1& Q0)const;
 public:
  int checkDensity()const;
  int applyDensityLimiter(double* min_density);
  double check_for_nonpos_density()const;
  //double get_minimum_density_at_positivity_points()const;
  void get_extreme_possible_states(
    int mcomponents, IntervalArray& Qinterval)const;
  CellPositivityLimiter(int i_in, int j_in, int moffset_in,
    dTensorBC4& q_in, const PositivityPointInfo& points_in):
   i(i_in), j(j_in), moffset(moffset_in), q(q_in), points(points_in)
  { }
};

class Gas05_CellPositivityLimiter: public CellPositivityLimiter
{
 public:
  double get_theta(const dTensor1& Q0,const dTensor2& dQ,
    const dTensor1& pos, int worst_offender)const;
  static double get_theta_for_point(int mp, const dTensor1& Q0,
    const dTensor2& dQ, double p0, double p1);
  int get_positivity_indicator_at_points(
    dTensor1& pos, const dTensor1& Q0,const dTensor2& dQ, double t,
    double epsilon)const;
  void reset_avg_thermal_energy(dTensor1& Q0, const double new_thermal_energy);
  int check()const;
  int apply();
  Gas05_CellPositivityLimiter(int i_in, int j_in, int moffset_in,
      dTensorBC4& q_in, const PositivityPointInfo& points) :
    CellPositivityLimiter(i_in,j_in,moffset_in,q_in,points) {}
};

class Offender;
class Gas10_CellPositivityLimiter: public CellPositivityLimiter
{
 private:
  void print_evals_and_pos_indicators_for_point(
    int mp,const dTensor1& Q0,const dTensor2& dQ,const dTensor2& pos)const;
  double get_theta_using_pos_indicators(
    Offender& worst_offender, dTensor2& pos,
    const dTensor1& Q0,const dTensor2& dQ,
    double epsilon)const;
  void shift_pressure(
    dTensor1& Q0, const double pshift);
 public:
  int check()const;
  int apply();
  Gas10_CellPositivityLimiter(int i_in, int j_in, int moffset_in,
    dTensorBC4& q_in, const PositivityPointInfo& points):
   CellPositivityLimiter(i_in,j_in,moffset_in,q_in,points)
  {}
};

namespace PositivityPointType{
  enum Enum {
    FLUX=1, // use Gauss-Lobatto points for flux limiters
    SOURCE=2, // use quadrature points for source limiters
    RIEMANN=3 // a subset of FLUX
  };
}

class PositivityPointInfo
{
 private:
  int numPoints;
  PositivityPointType::Enum pointType;
  const dTensor2* polynomValsAtPoints;
  const IntervalArray* polynomIntervals;
 private:
  void set_numPoints(int in){numPoints = in;}
  void set_polynomIntervals(const IntervalArray* in){polynomIntervals = in;}
  void set_polynomValsAtPoints(const dTensor2* in){polynomValsAtPoints = in;}
 public:
  PositivityPointType::Enum get_pointType()const{return pointType;}
  int get_numPoints() const {return numPoints;}
  const IntervalArray& get_polynomIntervals()const{return *polynomIntervals;}
  const dTensor2& get_polynomValsAtPoints()const{return *polynomValsAtPoints;}

  PositivityPointInfo(PositivityPointType::Enum pointType_in);
};

class PositivityLimiter
{
 private:
  PositivityPointInfo points;
  dTensorBC4& q;
 private:
  //void setPositivityPointInfo(PositivityPointType::Enum pointType);
  int ApplyInCell(int i,int j);
 public:
  int checkPositivity()const;
  int checkPositivityInCell(int i,int j) const;
  int apply();
  PositivityLimiter(
    dTensorBC4& q_in,
    PositivityPointType::Enum pointType):
   q(q_in),
   points(pointType)
  {};
};

//void ApplyPositivityLimiter(dTensorBC4& q);
bool gas10_verify_positivity_of_cell_averages(
  int moffset, const dTensorBC4& q, double* pos0min);
bool gas10_enforce_positivity_of_cell_averages(
  int moffset, dTensorBC4& q, double* pos0min);

#endif
