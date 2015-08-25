//#define MAX_DEBUG_LEVEL 3 //comment out to turn off verbose debug in this file
#include "debug.h"
#include <float.h>
#include "tensors.h"
#include "Limiters.h"
#include "Interval.h"
#include "Legendre2d.h"
#include "DogParams.h"
#include "FiveMomentParams.h"
#include "PlasmaParams.h"
#include "gas05.h"
#include "gas10.h"
#include "math.h"
#include "dog_math.h"
#include "constants.h"
#include "Polynomial.h"

// for the case of Interval checks using invariants is
// not clearly desirable and increases expense
//#define Interval_USE_INVARIANT_POS_IND
#define USE_INVARIANT_POS_IND

// *** begin positivity_limiters ***
using namespace PositivityReturnType;

const char* get_string(PositivityReturnType::Enum in)
{
  switch(in)
  {
    case NEG_PRESSURE:
      return "NEG_PRESSURE";
    case NEG_DENSITY:
      return "NEG_DENSITY";
    case ENFORCED_DENSITY:
      return "ENFORCED_DENSITY";
    case ENFORCED_PRESSURE:
      return "ENFORCED_PRESSURE";
    case ENFORCED_POSITIVITY:
      return "ENFORCED_POSITIVITY";
    case EVALUATIONS_CONFIRMED_POSITIVITY:
      return "EVALUATIONS_CONFIRMED_POSITIVITY";
    case INTERVALS_CONFIRMED_POSITIVITY:
      return "INTERVALS_CONFIRMED_POSITIVITY";
    case UNDEFINED:
      return "UNDEFINED";
  }
  invalid_value_error(in);
  return "error";
}

// =============== generic stuff ================

// we get a more precise estimate if we use interval
// arithmetic rather than infinity norm of the deviation.
//
void CellPositivityLimiter::get_extreme_possible_states(
  int mcomponents, IntervalArray& Qinterval)const
{
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const IntervalArray& phi_interval = points.get_polynomIntervals();
  // Legendre2d::get_phi_interval();
  const int mstart = 1+moffset;
  const int mend = mcomponents+moffset;
  for(int m=mstart;m<=mend;m++)
  {
    //int double q_max_deviation=0.;
    Interval q_dev(0.,0.);
    for(int k=2;k<=kmax;k++)
    {
      q_dev += phi_interval.get(k)*q.get(i,j,m,k);
      // q_max_deviation += fabs(q.get(i,j,m,k))*phi_max(k);
    }
    const double q_ave = q.get(i,j,m,1);
    Qinterval.fetch(m-moffset) = q_dev + q_ave;
  }
}

void CellPositivityLimiter::get_Q_at_point(dTensor1& Q, int mp)const
{
  const dTensor2& phi = points.get_polynomValsAtPoints();
  assert_ge(mp,1);
  assert_le(mp,phi.getsize(1));
  const int kmax=phi.getsize(2);
  for(int mQ=1;mQ<=Q.getsize();mQ++)
  {
    const int meq = mQ + moffset;
    double q_meq=q.get(i,j,meq,1); // initialize to average value
    for(int k=2;k<=kmax;k++)
    {
      // add deviations
      q_meq += q.get(i,j,meq,k)*phi.get(mp,k);
    }
    Q.set(mQ, q_meq);
  }
}

void CellPositivityLimiter::get_dQ_at_positivity_points(dTensor2& dQ)const
{
  const int mcomponents=dQ.getsize(2);
  //const dTensor2& phi = Legendre2d::get_phiAtPositivityPoints();
  const dTensor2& phi = points.get_polynomValsAtPoints();
  const int kmax=phi.getsize(2);
  for(int mp=1;mp<=dQ.getsize(1);mp++)
  for(int mQ=1;mQ<=mcomponents;mQ++)
  {
    const int meq = mQ + moffset;
    double dq=0.;
    for(int k=2;k<=kmax;k++) dq += q.get(i,j,meq,k)*phi.get(mp,k);
    dQ.set(mp,mQ, dq);
  }
}

inline void put_Q0_in_cell(int i,int j,int moffset,
  dTensorBC4& q, const dTensor1& Q0)
{
  const int mcomponents=Q0.getsize();
  for(int mQ=1;mQ<=mcomponents;mQ++)
  {
    const int meq = mQ + moffset;
    q.set(i,j,meq,1, Q0.get(mQ));
  }
}

inline void get_Q0_in_cell(int i,int j,int moffset,
  const dTensorBC4& q, dTensor1& Q0)
{
  const int mcomponents=Q0.getsize();
  for(int mQ=1;mQ<=mcomponents;mQ++)
  {
    const int meq = mQ + moffset;
    Q0.set(mQ, q.get(i,j,meq,1)); //*phi.get(mp,1), which equals 1.
  }
}

void CellPositivityLimiter::get_Q0_in_cell(dTensor1& Q0)const
{
  ::get_Q0_in_cell(i,j,moffset, q, Q0);
}

void get_Q_at_point(dTensor1& Q,
  int mp, const dTensor1& Q0, const dTensor2& dQ)
{
  for(int mQ=1;mQ<=dQ.getsize(2);mQ++)
  {
    Q.set(mQ, Q0.get(mQ)+dQ.get(mp,mQ));
  }
}

void get_Q_for_t_at_point(dTensor1& Q,
  int mp, const dTensor1& Q0, const dTensor2& dQ, double t)
{
  for(int mQ=1;mQ<=dQ.getsize(2);mQ++)
  {
    Q.set(mQ, Q0.get(mQ)+t*dQ.get(mp,mQ));
  }
}

void rescale_dq_in_cell(int i,int j,int moffset,int mcomponents,
  dTensorBC4& q, double theta)
{
  const int kmax = q.getsize(4);
  const int mstart = moffset+1;
  const int mend = moffset+mcomponents;
  for(int meq=mstart;meq<=mend;meq++)
  for(int k=2;k<=kmax;k++)
  {
    q.vfetch(q.getidx(i,j,meq,k)) *= theta;
  }
}

// =============== density stuff ================
double CellPositivityLimiter::get_minimum_density_at_points(
  const dTensor2& phi)const
{
  const int kmax = dogParams.get_kmax();
  const int meq = moffset+1; // assume density is first component
  double rho_min=DBL_MAX;
  for(int mp=1;mp<=phi.getsize(1);mp++)
  {
    double rho=0;
    // evaluate at point
    for(int k=1;k<=kmax;k++) rho += q.get(i,j,meq,k)*phi.get(mp,k);
    if(rho<rho_min)
    {
      rho_min = rho;
      //mpoint=mp;
    }
  }
  return rho_min;
}

//double CellPositivityLimiter::get_minimum_density_at_coarse_points()const
//{
//  return get_minimum_density_at_points(Legendre2d::get_phiCoarse());
//}

//double CellPositivityLimiter::get_minimum_density_at_positivity_points()const
//{
//  // iterate over set of points where we enforce positivity
//  //const dTensor2& phi = Legendre2d::get_phiAtPositivityPoints();
//  return get_minimum_density_at_points(points.get_polynomValsAtPoints());
//}

// if density falls below limit return minimum density among
// positivity points; else return a lower bound on density
// at those points
double CellPositivityLimiter::check_for_nonpos_density()const
{
  //double rho_min[2];
  IntervalArray rho(2);
  get_extreme_possible_states(1,rho);
  const double rho_min = rho.get(1).min();
  if(rho_min >= LimiterParams::get_rho_epsilon())
    return rho_min;

  // this point in the code is rarely reached
  dprintf3("for (%d,%d,%d) rho_min = %24.16e",i,j,moffset,rho_min);

  // get minimum density at positivity points
  double min_density
    = get_minimum_density_at_points(points.get_polynomValsAtPoints());
  // lest the density become zero at a flux point
  if(points.get_pointType()==PositivityPointType::FLUX)
  {
    // get minimum density at coarse points
    min_density = Min(min_density,
      get_minimum_density_at_points(Legendre2d::get_phiCoarse()));
  }
  return min_density;
}

int CellPositivityLimiter::checkDensity()const
{
  const double rho_epsilon = LimiterParams::get_rho_epsilon();

  const double rho_ave = q.get(i,j,moffset+1,1);
  if(rho_ave < rho_epsilon)
  {
    dprintf1("cell average density violates positivity: %24.16e",rho_ave);
    return NEG_DENSITY;
  }

  IntervalArray rho(2);
  get_extreme_possible_states(1,rho);
  if(rho.get(1).min() >= rho_epsilon)
    return INTERVALS_CONFIRMED_POSITIVITY;

  double min_density = check_for_nonpos_density();
  if(min_density >= rho_epsilon)
    return EVALUATIONS_CONFIRMED_POSITIVITY;

  dprint1(min_density);
  return NEG_DENSITY;
}

// returns lower bound on density
int CellPositivityLimiter::applyDensityLimiter(double* min_density)
{
  using namespace LimiterParams;
  *min_density = check_for_nonpos_density();
  if(*min_density >= get_rho_epsilon())
    return EVALUATIONS_CONFIRMED_POSITIVITY;

  // this point in the code is only reached
  // on the very rare occasion that the density
  // is not positive at a positivity point
  //
  dprintf2("for (%d,%d,%d) *min_density = %24.16e",
    i,j,moffset,*min_density);

  // rescale legendre coefficients to make density positive
  // at all the points we tested.
  int meq = moffset+1;
  double rho_0 = q.get(i,j,meq,1);
  //
  // ensure that average density is positive
  //double min_density = ensurePositiveAverageDensity();
  if(rho_0 < 0)
  {
    eprintf("positivity violation: rho_0 = %24.16e",rho_0);
    //timeStepException.Throw(NEGATIVE_DENSITY);
  }
  else if(rho_0 < LimiterParams::get_rho_epsilon())
  {
    // shift solution by a constant so that rho_0 becomes epsilon
    const double rho_old = rho_0;
    rho_0 = LimiterParams::get_rho_epsilon();
    double shift = rho_0-rho_old;
    *min_density += shift;
    q.set(i,j,1+moffset,1, rho_0);
  }
  assert_ge(rho_0,get_rho_epsilon());
  double theta = (rho_0-get_rho_epsilon())/(rho_0-*min_density);
  assert_le(theta, 1.);

  // limit deviations of q
  for( int k=2; k <= dogParams.get_kmax(); k++ )
    q.set(i,j,meq,k, q.get(i,j,meq,k) * theta);

  // perform check
  if(true)
  {
    *min_density = check_for_nonpos_density();
    assert_ge(*min_density,get_rho_epsilon()*0.5);
  }
  return ENFORCED_DENSITY;
}

// =============== gas05 stuff ================

// return the machine pressure positivity indicator, which is the
// thermal energy minus the epsilon for the pressure indicator:
//
//   p0 = thermal_energy - epsilon
//
// want p/(gamma-1) >= eps, i.e.
// (N-0.5*M^2/rho) >= eps, i.e.
// (N-eps-0.5*M^2/rho) >= 0, i.e.
// (rho*(N-eps)-0.5*M^2) >= 0,
// which is a quadratic inequality in t.
//
double gas05_get_pos_indicator_for_cell_avg_state(
  const dTensor1& Q0,
  const double epsilon)
{
  using namespace FiveMomentComponentID;
  const double rho = Q0.get(_rho);
  const double M1  = Q0.get(_M1);
  const double M2  = Q0.get(_M2);
  const double M3  = Q0.get(_M3);
  const double N   = Q0.get(_N);
  return rho*(N-epsilon)-0.5*(M1*M1+M2*M2+M3*M3);
}

// this should only be called to increase the thermal energy
// to a machine minimum
void Gas05_CellPositivityLimiter::reset_avg_thermal_energy(
  dTensor1& Q0,
  const double new_thermal_energy)
{
  using namespace FiveMomentComponentID;
  const double rho = Q0.get(_rho);
  const double M1  = Q0.get(_M1);
  const double M2  = Q0.get(_M2);
  const double M3  = Q0.get(_M3);
  const double N   = Q0.get(_N);
  const double KE = 0.5*(M1*M1+M2*M2+M3*M3)/rho;
  const double thermal_energy = N-KE;
  if(thermal_energy<0.)
  {
    eprintf("positivity violation: thermal_energy=%24.16e\n");
    //timeStepException.Throw(NEGATIVE_PRESSURE);
  }
  assert_le(thermal_energy, new_thermal_energy);
  const double new_energy = new_thermal_energy + KE;
  Q0.set(_N,new_energy);
  q.set(i,j,moffset+_N,1, new_energy);
}


inline double get_gas05_positivity_indicator_at_point(
  const dTensor1& Q0, const dTensor2& dQ, int mp, double t, double epsilon)
{
  using namespace FiveMomentComponentID;
  dTensor1 Q(5);
  get_Q_for_t_at_point(Q, mp,Q0,dQ,t);
  const double rho = Q.get(_rho);
  const double M1  = Q.get(_M1);
  const double M2  = Q.get(_M2);
  const double M3  = Q.get(_M3);
  const double N   = Q.get(_N);

  return rho*(N-epsilon)-0.5*(M1*M1+M2*M2+M3*M3);
}

// get the pressure at all the points; return the
// index of a worst offender (a point with minimum pressure)
// or 0 if there is no offender
int Gas05_CellPositivityLimiter::get_positivity_indicator_at_points(
  dTensor1& pos, const dTensor1& Q0,const dTensor2& dQ, double t,
  double epsilon)const
{
  int min_point = 0;
  // iterate over set of points where we enforce positivity
  const int kmax = dogParams.get_kmax();
  //const dTensor2& phi = Legendre2d::get_phiAtPositivityPoints();
  const dTensor2& phi = points.get_polynomValsAtPoints();
  double pos_min=DBL_MAX;
  for(int mp=1;mp<=phi.getsize(1);mp++)
  {
    double pos_indicator = get_gas05_positivity_indicator_at_point(
      Q0,dQ,mp,t,epsilon);
    pos.set(mp,pos_indicator);
    if(pos_indicator<pos_min)
    {
      pos_min = pos_indicator;
      min_point=mp;
    }
  }
  if(pos_min >= 0.) return 0;
  return min_point;
}

 // return smallest positive theta such that
 // positivity_indicator(theta) = 0.
 // else return 1.;
double Gas05_CellPositivityLimiter::get_theta_for_point(int mp,
  const dTensor1& Q0,
  const dTensor2& dQ, double p0, double p1)
{
   double theta = 1.;
   const double epsilon = LimiterParams::get_pressure_epsilon();
   double ph = get_gas05_positivity_indicator_at_point(
     Q0, dQ, mp, 0.5, epsilon);

   // determine coefficients of quadratic p(t) = at^2+b*t+c
   // given that p(0) = p0, p(0.5) = ph, and p(1) = p1
   //const double ph = p_half.get(mp);
   const double& c = p0;
   const double b = 4*ph-p1-3*p0;
   const double a = p1-b-p0;
   const double bh = b/2.;
   const double discriminant = bh*bh-a*c;
   if(discriminant<0.) return 1.; // irreducible quadratic
   if(fabs(a)<EPSILON) // linear case
   {
     if(fabs(bh)<EPSILON) // constant case
     {
       if(fabs(c)<EPSILON) // zero case (should never get here anyway)
       {
         return 0.; // all deviations should be zeroed.
       }
       else
       {
         return 1.; // constant case
       }
     }
     const double theta = -c/b; // linear case
   }
   else // parabola with roots
   {
     const double sq_disc = sqrt(discriminant);
     double t1 = (-bh - sq_disc)/a;
     double t2 = (-bh + sq_disc)/a;
     double tmin,tmax;
     if(t1<t2) {
       tmin = t1;
       tmax = t2;
     }
     else
     {
       tmin = t2;
       tmax = t1;
     }
     if(tmin>=0.) theta=tmin;
     else if(tmax>=0) theta=tmax;
   }
   return theta;
}

double Gas05_CellPositivityLimiter::get_theta(
  const dTensor1& Q0, 
  const dTensor2& dQ,
  const dTensor1& pos,
  int worst_offender)const
{
  const double epsilon = LimiterParams::get_pressure_epsilon();
  double p0 = gas05_get_pos_indicator_for_cell_avg_state(Q0,epsilon);
  //dTensor1 p_half(numPoints); // p(0.5)
  //get_positivity_indicator_at_points(p_half, Q0, dQ, 0.5);

  // for each point solve p(t)=0 for 0<t<1
  // and take the minimum of all such t.
  double theta = 1;
  int true_worst_offender = 0;
  //const int numPoints = Legendre2d::get_numPositivityPoints();
  const int numPoints = points.get_numPoints();
  for(int mp=1;mp<=numPoints;mp++)
  {
    const double p1 = pos.get(mp);
    double point_theta = get_theta_for_point(mp, Q0, dQ, p0, p1);
    if(point_theta > 0. && point_theta < theta)
    {
      theta = point_theta;
      true_worst_offender = mp;
    }
  }
  // maintain statistic on whether true_worst_offender==worst_offender ?

  // check theta
  assert_le(theta,1.);
  assert_le(0.,theta);
}

// we want to enforce that thermal energy is greater than
// or equal to epsilon
Interval gas05_get_pos_ind_interval(const IntervalArray& Q, double epsilon)
{
  using namespace FiveMomentComponentID;
  const Interval& rho = Q.get(_rho);
  const Interval& M1 = Q.get(_M1);
  const Interval& M2 = Q.get(_M2);
  const Interval& M3 = Q.get(_M3);

  const Interval rho_KE = 0.5*(M1.square()+M2.square()+M3.square());
  const Interval pos = rho*(Q.get(_N)-epsilon) - rho_KE;
  return pos;
}

int Gas05_CellPositivityLimiter::check()const
{
  using namespace PositivityReturnType;
  int ret = MAX;

  ret = iMin(ret, checkDensity());
  // not much point in checking the pressure if the density is negative
  if(ret == NEG_DENSITY) return ret;

  // remainder checks positivity of pressure
  
  // quick check to try to verify that pressure is positive
  //
  IntervalArray Qinterval(6);
  get_extreme_possible_states(5,Qinterval);
  // since positivity limiting has been applied to the density we in fact have:
  Interval pressure = gas05_get_pos_ind_interval(Qinterval,
    LimiterParams::get_pressure_epsilon());
  // verified positive?
  if(pressure.min() >= LimiterParams::get_pressure_epsilon())
  {
    ret = iMin(ret, INTERVALS_CONFIRMED_POSITIVITY);
    return ret;
  }

  // do_we_have_positivity_at_all_the_positivity_points?
  // (intervals were not enough, so check state at positivity points.)
  //
  const int numPoints = points.get_numPoints();
  dTensor1 Q0(5);
  get_Q0_in_cell(Q0);
  double p0 = gas05_get_pos_indicator_for_cell_avg_state(
    Q0,LimiterParams::get_avg_pressure_epsilon());
  if(p0<=0.)
  {
    dprintf1("cell average violates positivity of pressure: %24.16e", p0);
    ret = iMin(ret, NEG_PRESSURE);
    return ret;
  }
  //
  dTensor2 dQ(numPoints,5);
  get_dQ_at_positivity_points(dQ);
  // check pressure positivity at all points
  dTensor1 pos(numPoints);
  const double epsilon = LimiterParams::get_pressure_epsilon();
  int worst_offender = get_positivity_indicator_at_points(pos,Q0,dQ,1.,epsilon);
  if(!worst_offender)
  {
    ret = iMin(ret, EVALUATIONS_CONFIRMED_POSITIVITY);
    return ret;
  }
  
  dprintf1("violation of positivity of pressure");
  ret = iMin(ret, NEG_PRESSURE);
  return ret;
}

int Gas05_CellPositivityLimiter::apply()
{
  int ret = PositivityReturnType::MAX;
  // enforce positivity of density
  double min_density;
  ret = iMin(ret,applyDensityLimiter(&min_density));

  // remainder enforces positivity of pressure
  
  const double pressure_epsilon = LimiterParams::get_pressure_epsilon();
  // the value that the average thermal energy is held above
  const double avg_eps = LimiterParams::get_avg_pressure_epsilon();

  // quick check to try to verify that pressure is positive
  //
  IntervalArray Qinterval(6);
  get_extreme_possible_states(5,Qinterval);
  // since positivity limiting has been applied to the density we in fact have:
  Qinterval.fetch(1).set_min(min_density);
  Interval pressure = gas05_get_pos_ind_interval(Qinterval,avg_eps);
  // verified positive?
  if(pressure.min() >= pressure_epsilon)
  {
    ret = iMin(ret, INTERVALS_CONFIRMED_POSITIVITY);
    return ret;
  }

  // do_we_have_positivity_at_all_the_positivity_points?
  // (intervals were not enough, so check state at positivity points.)
  //
  const int numPoints = points.get_numPoints();
  //Legendre2d::get_numPositivityPoints();
  dTensor1 Q0(5);
  get_Q0_in_cell(Q0);
  // ensure that positivity indicator of cell average is strictly positive
  double p0 = gas05_get_pos_indicator_for_cell_avg_state(Q0,avg_eps);
  if(p0<=0.)
  {
    double p00 = gas05_get_pos_indicator_for_cell_avg_state(Q0,0.);
    assert_gt(p00,0.);

    // increase the thermal energy to exactly avg_eps,
    reset_avg_thermal_energy(Q0, avg_eps);
    // in this case one is tempted to set all higher polynomial
    // moments to zero and return, but this is not continuous
    // so instead we proceed.
    // (but there is no requirement to be continuous on the
    // order of machine epsilon...)
    //if(avg_eps <= pressure_epsilon+EPSILON)
    //{
    //  rescale_dq_in_cell(i,j,moffset,5,q,0.); // assign q := q0
    //  return ENFORCED_POSITIVITY;
    //}
  }
  //
  dTensor2 dQ(numPoints,5);
  get_dQ_at_positivity_points(dQ);
  // check pressure positivity at all points
  dTensor1 pos(numPoints);
  int worst_offender = get_positivity_indicator_at_points(
    pos,Q0,dQ,1.,pressure_epsilon);
  if(!worst_offender)
  {
    ret = iMin(ret, EVALUATIONS_CONFIRMED_POSITIVITY);
    return ret;
  }
  
  // enforce positivity at all points
  //
  // Let q = q0 + t*(q-q0), where q0 is cell average.
  // at each (positivity) point the pressure p is a quadratic polynomial in t.
  // Find theta = the maximum t between 0 and 1
  // such that for all t<theta p(t)>=epsilon for all positivity points.
  const double theta = get_theta(Q0,dQ,pos,worst_offender);
  rescale_dq_in_cell(i,j,moffset,5,q,theta); // assign q := q0+theta(q-q0)

  // confirm positivity at all points
  // (comment this out once the code is working)
  #if 1
  get_dQ_at_positivity_points(dQ);
  double p_min = get_positivity_indicator_at_points(
    pos,Q0,dQ,1.,0.5*pressure_epsilon);
  assert_ge(p_min,0.);
  #endif
  ret = iMin(ret, ENFORCED_PRESSURE);
  return ret;
}

// =============== gas10 stuff ================

// positivity indicator index
class PII
{
  int point_idx;
  int poly_idx;
 public:
  PII():point_idx(0),poly_idx(0){};
  int get_point_idx()const{return point_idx;}
  int get_poly_idx ()const{return poly_idx ;}
  void set_poly_idx(int poly_idx_){poly_idx=poly_idx_;}
  void set(int point_idx_, int poly_idx_)
  {
    point_idx=point_idx_;
    poly_idx=poly_idx_;
  }
};

class Offender : public PII
{
  double severity;
  double value;
 private:
  // disabled:
  void set(int point_idx_, int poly_idx_);
  void set_poly_idx(int poly_idx_);
 public:
  Offender():severity(0.),value(0.){}
  bool not_set()const{return severity<=0.;}
  bool is_set()const{return severity>0.;}
  double get_severity(){return severity;}
  void set(int point_idx_, int poly_idx_, double severity_, double value_){
    PII::set(point_idx_, poly_idx_);
    severity=severity_;
    value = value_;
  }
  //void set_poly_idx(int poly_idx_,double severity_,double value_){
  //  PII::set_poly_idx(poly_idx_);
  //  severity=severity_;
  //  value=value_;
  //}
  void print()
  {
    printf("%d.%d: severity=%24.16e, value=%24.16e",
      get_point_idx(),get_poly_idx(),severity,value);
  }
};

void gas10_rescale_pressure_anisotropy(dTensor1& prim,
  double rescaling_factor)
{
  using namespace TenMomentComponentID;
  const double P11 = prim.get(_P11);
  const double P22 = prim.get(_P22);
  const double P33 = prim.get(_P33);
  const double p = (1./3.)*(P11+P22+P33);
  const double dP11 = P11-p;
  const double dP22 = P22-p;
  const double dP33 = P33-p;

  prim.set(_P11, p+dP11*rescaling_factor);
  prim.set(_P22, p+dP22*rescaling_factor);
  prim.set(_P33, p+dP33*rescaling_factor);

  prim.fetch(_P12) *= rescaling_factor;
  prim.fetch(_P13) *= rescaling_factor;
  prim.fetch(_P23) *= rescaling_factor;
}

double gas10_get_pressure_evals_from_prim(const dTensor1& prim, double* evals)
{
  using namespace TenMomentComponentID;
  const double P11 = prim.get(_P11);
  const double P12 = prim.get(_P12);
  const double P13 = prim.get(_P13);
  const double P22 = prim.get(_P22);
  const double P23 = prim.get(_P23);
  const double P33 = prim.get(_P33);

  gas10_get_pressure_evals(evals, P11,P12,P13,P22,P23,P33);
  return (1./3.)*(P11+P22+P33);
}

void gas10_get_pressure_evals(const dTensor1& Q0, double* evals)
{
  using namespace TenMomentComponentID;
  const double rho  = Q0.get(_rho);
  const double M1   = Q0.get(_M1);
  const double M2   = Q0.get(_M2);
  const double M3   = Q0.get(_M3);
  const double Ne11 = Q0.get(_N11);
  const double Ne12 = Q0.get(_N12);
  const double Ne13 = Q0.get(_N13);
  const double Ne22 = Q0.get(_N22);
  const double Ne23 = Q0.get(_N23);
  const double Ne33 = Q0.get(_N33);

  const double rho_inv = 1./rho;
  const double P11 = Ne11 - M1*M1*rho_inv;
  const double P12 = Ne12 - M1*M2*rho_inv;
  const double P13 = Ne13 - M1*M3*rho_inv;
  const double P22 = Ne22 - M2*M2*rho_inv;
  const double P23 = Ne23 - M2*M3*rho_inv;
  const double P33 = Ne33 - M3*M3*rho_inv;

  gas10_get_pressure_evals(evals, P11,P12,P13,P22,P23,P33);
}

void Gas10_CellPositivityLimiter::shift_pressure(
  dTensor1& Q0, const double pshift)
{
  using namespace TenMomentComponentID;
  // shift the pressure tensor of the cell average
  const double N11 = Q0.get(_N11);
  const double N22 = Q0.get(_N22);
  const double N33 = Q0.get(_N33);
  const double N11_new = N11+pshift;
  const double N22_new = N22+pshift;
  const double N33_new = N33+pshift;
  Q0.set(_N11,N11_new);
  Q0.set(_N22,N22_new);
  Q0.set(_N33,N33_new);
  q.set(i,j,moffset+_N11,1,  N11_new);
  q.set(i,j,moffset+_N22,1,  N22_new);
  q.set(i,j,moffset+_N33,1,  N33_new);
}

// first use Gershgorin (hoping for a faster estimate that
// applies to most cases) (should do statistics to see how
// frequently this test returns true)
inline bool gas10_guaranteed_positivity_using_Gershgorin(
  const Interval& Pos11,
  const Interval& Pos12,
  const Interval& Pos13,
  const Interval& Pos22,
  const Interval& Pos23,
  const Interval& Pos33)
{
  double Pos12abs = Max(fabs(Pos12.min()),fabs(Pos12.max()));
  double Pos13abs = Max(fabs(Pos13.min()),fabs(Pos13.max()));
  double Pos23abs = Max(fabs(Pos23.min()),fabs(Pos23.max()));
  double pmin = Pos11.min() - Pos12abs - Pos13abs;
  pmin = Min(pmin,Pos22.min() - Pos12abs - Pos23abs);
  pmin = Min(pmin,Pos33.min() - Pos13abs - Pos23abs);
  if(pmin>=0.)
  {
    return true;
  }
  return false;
}

// ensure that thermal energy tensor exceeds epsilon
// (that is, TE - eps*identity_tensor is positive definite)
bool gas10_guaranteed_positivity_using_intervals(
  const IntervalArray& Q,
  double epsilon)
{
  using namespace TenMomentComponentID;
  const Interval& rho = Q.get(_rho);
  const Interval& M1 = Q.get(_M1);
  const Interval& M2 = Q.get(_M2);
  const Interval& M3 = Q.get(_M3);

  const Interval Ne11 = Q.get(_N11) - epsilon;
  const Interval Ne12 = Q.get(_N12);
  const Interval Ne13 = Q.get(_N13);
  const Interval Ne22 = Q.get(_N22) - epsilon;
  const Interval Ne23 = Q.get(_N23);
  const Interval Ne33 = Q.get(_N33) - epsilon;

  const Interval Pos11 = rho*Ne11 - M1.square();
  const Interval Pos12 = rho*Ne12 - M1*M2;
  const Interval Pos13 = rho*Ne13 - M1*M3;
  const Interval Pos22 = rho*Ne22 - M2.square();
  const Interval Pos23 = rho*Ne23 - M2*M3;
  const Interval Pos33 = rho*Ne33 - M3.square();

  // should take this out unless it returns true more than
  // say 70% of the time
  //
  #if 0 // never seems to return true so why bother
  if(gas10_guaranteed_positivity_using_Gershgorin(
     Pos11, Pos12, Pos13, Pos22, Pos23, Pos33))
  {
    return true;
  }
  #endif

  // return true only if the pressure is guaranteed positive-definite
  //
#ifndef Interval_USE_INVARIANT_POS_IND
  if(!(Pos11 >= 0.))
#else /* Interval_USE_INVARIANT_POS_IND */
  //if(!(Pos11 >= 0.))
  Interval trPos = Pos11+Pos22+Pos33;
  if(!(trPos >= 0.))
#endif /* Interval_USE_INVARIANT_POS_IND */
  {
    //dprint3(Pos11.min());
    return false;
  }
  else
  {
    Interval Adj33 = Pos11*Pos22-Pos12*Pos12;
#ifndef Interval_USE_INVARIANT_POS_IND
    if(!(Adj33>=0.))
#else /* Interval_USE_INVARIANT_POS_IND */
    // extra computation to use matrix invariants (needed for stability?)
    Interval Adj11 = Pos22*Pos33 - Pos23*Pos23;
    Interval Adj22 = Pos33*Pos11 - Pos13*Pos13;
    Interval trAdj = Adj11+Adj22+Adj33;
    //if(!(Adj33>=0.))
    if(!(trAdj>=0.))
#endif /* Interval_USE_INVARIANT_POS_IND */
    {
      //dprint3(Adj33.min());
      return false;
    }
    else
    {
      // compute components of adjugate needed for determinant.
      Interval Adj13 = Pos12*Pos23 - Pos22*Pos13;
      Interval Adj23 = Pos12*Pos13 - Pos11*Pos23;
      Interval detP = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;
      if(!(detP>=0.))
      {
        //dprint3(detP.min());
        return false;
      }
    }
  }

  return true;
}

// Use q to denote conserved variables.
// Let q(t) = q0 + t*(q1-q0),
// where q0 is cell average,
// where q1 is value (at some positivity point), and
// t = real number between 0 and 1.
//
// For brevity of notation adopt the following matrix conventions:
// (1) A > B if (A-B) > 0, i.e. if (A-B) is strictly positive definite and
//     A >= B if (A-B) >= 0, i.e. if (A-B) is non-strictly positive definite;
// (2) A + s = A + s*I.
//
// We wish to enforce the positivity condition
// Pr >= eps, where eps is a "positivity epsilon" and Pr is Pressure
// (which will be significantly greater than machine epsilon).
// That is, we want Pr - eps >= 0.  But Pr = E-M*M/rho.
// We want rho*(Pr-eps) >= 0.  That is, rho*(E-eps)-M*M >= 0.
// Let P_eps = rho*(E-eps)-M*M.
//
// Let P = P_0 = rho*E-M*M = P_eps + rho*eps.
// Write reps = rho*eps.  Then
// det P = det(P_eps) + reps*tr(cof(P_eps)) + reps^2*tr(P_eps) + reps^3
// = Le1*Le2*Le3 + reps*(Le1*Le3+Le1*Le2+Le2*Le3)+reps^2(Le1+Le2+Le3) + reps^3,
// where Le1 <= Le2 <= Le3 are the eigenvalues of P_eps.
// 
// We can use a root-finder to enforce that
// (1) Le1*Le2*Le3 = det(P_eps) >= Order(eps_mach),
// (2) (Le1*Le3+Le1*Le2+Le2*Le3) = tr(cof(P_eps)) >= Order(eps_mach), and
// (3) (Le1+Le2+Le3) = tr(P_eps) >= Order(eps_mach).
//
// Let L1 <= L2 <= L3 be the eigenvalues of P.
// Note that L1 = Le1 + rho*eps, etc.
// We want that Le1 >= Order(eps_mach).
// Then we would have that
// L1 >= rho*eps + Order(eps_mach).
//
// Want to ensure that the condition enforced by the root-finder
// makes it impossible to have Le1 < -Order(eps_mach).
//
//
// For all directions n such that n.n = 1 we need
// n.P.n >= eps, i.e.
// n.(E-M*M/rho).n >= eps, i.e.
// n.(E-M*M/rho).n >= n.eps*I.n, i.e.
// n.(E-eps*I-M*M/rho).n >= 0, i.e.
// n.(rho*(E-eps*I)-M*M).n >= 0, i.e.
//   (rho*(E-eps*I)-M*M) is positive definite, i.e.
// det(Pos^(n)) >= 0 for all principle square matrices Pos^(n)
//   of Pos(t) := (rho*(E-eps*I)-M*M);
// 
// which gives a quadratic, a quartic, and a sixth-order polynomial
// inequality in t.
//
// For numerical stability our positivity indicators should be
// invariants: the trace, the sum of the determinants of
// minor matrices of dimension two, i.e. the trace of the cofactor
// matrix, i.e. the sum of the products of pairs of eigenvalues,
// and the determinant. Descarte's rule of signs says that if
// these are all positive then the eigenvalues themselves must be
// positive.
// This is because guaranteeing that the product of eigenvalues
// is small does not guarantee that the eigenvalues themselves
// are small.
//
// Let t_{max} = argmin_{t in [0,1]} Pos(t) is not strictly positive definite.
// Observe that t_{max} = argmin_{t in [0,1]} det(Pos(t)) = 0.
// So it is enough to find the smallest t such that det(Pos(t)) = 0.
// 
bool get_gas10_positivity_indicators_for_cell_avg_state(
  const dTensor1& Q0, double* pos0, double epsilon)
{
  using namespace TenMomentComponentID;
  using namespace LimiterParams;
  const double rho  = Q0.get(_rho);
  const double M1   = Q0.get(_M1);
  const double M2   = Q0.get(_M2);
  const double M3   = Q0.get(_M3);
  const double Ne11 = Q0.get(_N11)-epsilon;
  const double Ne12 = Q0.get(_N12);
  const double Ne13 = Q0.get(_N13);
  const double Ne22 = Q0.get(_N22)-epsilon;
  const double Ne23 = Q0.get(_N23);
  const double Ne33 = Q0.get(_N33)-epsilon;

  const double Pos11 = rho*Ne11 - M1*M1;
  const double Pos12 = rho*Ne12 - M1*M2;
  const double Pos13 = rho*Ne13 - M1*M3;
  const double Pos22 = rho*Ne22 - M2*M2;
  const double Pos23 = rho*Ne23 - M2*M3;
  const double Pos33 = rho*Ne33 - M3*M3;

  // compute components of adjugate needed for determinant.
  const double Adj13 = Pos12*Pos23 - Pos22*Pos13;
  const double Adj23 = Pos12*Pos13 - Pos11*Pos23;
  const double Adj33 = Pos11*Pos22 - Pos12*Pos12;

#ifdef USE_INVARIANT_POS_IND
  // extra computation so as to use matrix invariants
  //
  const double Adj11 = Pos22*Pos33 - Pos23*Pos23;
  const double Adj22 = Pos33*Pos11 - Pos13*Pos13;
  const double trAdj = Adj11+Adj22+Adj33;
  const double trPos = Pos11+Pos22+Pos33;

#endif /* USE_INVARIANT_POS_IND */
  const double detPos = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;

  pos0[0]=rho-get_rho_epsilon();
  // positive-definite means that all three of these are positive.
#ifndef USE_INVARIANT_POS_IND
  pos0[1]=Pos11;
  pos0[2]=Adj33;
#else /* USE_INVARIANT_POS_IND */
  pos0[1]=trPos;
  pos0[2]=trAdj;
#endif /* USE_INVARIANT_POS_IND */
  pos0[3]=detPos;

  if(pos0[0]<0.||pos0[1]<0.||pos0[2]<0.||pos0[3]<0.) return false;

  return true;
}

double get_gas10_positivity_indicator_at_point(int mp, int mi,
  const dTensor1& Q0, const dTensor2& dQ, double t, double* pos=NULL)
{
  using namespace TenMomentComponentID;
  const double epsilon = LimiterParams::get_pressure_epsilon();
  dTensor1 Q(10);
  get_Q_for_t_at_point(Q, mp,Q0,dQ,t);
  const double rho = Q.get(_rho);
  const double M1  = Q.get(_M1);
  const double Ne11 = Q.get(_N11)-epsilon;
#ifdef USE_INVARIANT_POS_IND
  const double Ne22 = Q.get(_N22)-epsilon;
  const double Ne33 = Q.get(_N33)-epsilon;
  const double M2  = Q.get(_M2);
  const double M3  = Q.get(_M3);
#endif /* USE_INVARIANT_POS_IND */
  const double Pos11 = rho*Ne11 - M1*M1;
#ifdef USE_INVARIANT_POS_IND
  const double Pos22 = rho*Ne22 - M2*M2;
  const double Pos33 = rho*Ne33 - M3*M3;
#endif /* USE_INVARIANT_POS_IND */
  //if(pos) pos[1] = Pos11;
#ifndef USE_INVARIANT_POS_IND
  if(mi<=1) return Pos11;
#else /* USE_INVARIANT_POS_IND */
  //if(mi<=1) return Pos11;
  if(mi<=1) return Pos11+Pos22+Pos33;
#endif /* USE_INVARIANT_POS_IND */

#ifndef USE_INVARIANT_POS_IND
  const double M2  = Q.get(_M2);
#endif /* ! USE_INVARIANT_POS_IND */
  const double Ne12 = Q.get(_N12);
#ifndef USE_INVARIANT_POS_IND
  const double Ne22 = Q.get(_N22)-epsilon;
#endif /* ! USE_INVARIANT_POS_IND */
  const double Pos12 = rho*Ne12 - M1*M2;
#ifndef USE_INVARIANT_POS_IND
  const double Pos22 = rho*Ne22 - M2*M2;
#endif /* ! USE_INVARIANT_POS_IND */
  const double Adj33 = Pos11*Pos22 - Pos12*Pos12;
#ifndef USE_INVARIANT_POS_IND
  //if(pos) pos[2] = Adj33;
  if(mi<=2) return Adj33;

  const double M3  = Q.get(_M3);
#endif /* ! USE_INVARIANT_POS_IND */
  const double Ne13 = Q.get(_N13);
  const double Ne23 = Q.get(_N23);
#ifndef USE_INVARIANT_POS_IND
  const double Ne33 = Q.get(_N33)-epsilon;

#endif /* ! USE_INVARIANT_POS_IND */
  const double Pos13 = rho*Ne13 - M1*M3;
  const double Pos23 = rho*Ne23 - M2*M3;
#ifndef USE_INVARIANT_POS_IND
  const double Pos33 = rho*Ne33 - M3*M3;
#else /* USE_INVARIANT_POS_IND */
  // extra computation so as to use matrix invariants
  const double Adj11 = Pos22*Pos33 - Pos23*Pos23;
  const double Adj22 = Pos33*Pos11 - Pos13*Pos13;
  //if(pos) pos[2] = Adj33;
  //if(mi<=2) return Adj33;
  if(mi<=2) return Adj11+Adj22+Adj33;
#endif /* USE_INVARIANT_POS_IND */

  const double Adj13 = Pos12*Pos23 - Pos22*Pos13;
  const double Adj23 = Pos12*Pos13 - Pos11*Pos23;

  const double detPos = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;
  //if(pos) pos[3] = detPos;
  return detPos;
}

// return 0 if okay, else return the index of the
// point that is worst-offending in some sense.
// (my hope is that most often there will be only one
// offender or one flagrant offender).
//
Offender get_gas10_positivity_indicators_at_points(dTensor2& pos,
  const dTensor1& Q0, const dTensor2& dQ, double t,double epsilon)
{
  using namespace TenMomentComponentID;
  Offender retval;

  // need positivity indicators of cell average to calculate worst offender
  //
  // confirm positivity indicators at cell average state
  double pos0[4];
  bool okay = get_gas10_positivity_indicators_for_cell_avg_state(
    Q0,pos0,epsilon);
  assert(okay); // this should have been taken care earlier

  // worst offender has maximal offense value
  //
  double max_offense_severity=-DBL_MAX;
  for(int mp=1;mp<=dQ.getsize(1);mp++)
  {
    dTensor1 Q(10);
    get_Q_for_t_at_point(Q, mp,Q0,dQ,t);
    const double rho = Q.get(_rho);
    const double M1  = Q.get(_M1);
    const double M2  = Q.get(_M2);
    const double M3  = Q.get(_M3);
    const double Ne11 = Q.get(_N11)-epsilon;
    const double Ne12 = Q.get(_N12);
    const double Ne13 = Q.get(_N13);
    const double Ne22 = Q.get(_N22)-epsilon;
    const double Ne23 = Q.get(_N23);
    const double Ne33 = Q.get(_N33)-epsilon;

    const double Pos11 = rho*Ne11 - M1*M1;
    const double Pos12 = rho*Ne12 - M1*M2;
    const double Pos13 = rho*Ne13 - M1*M3;
    const double Pos22 = rho*Ne22 - M2*M2;
    const double Pos23 = rho*Ne23 - M2*M3;
    const double Pos33 = rho*Ne33 - M3*M3;

    const double Adj13 = Pos12*Pos23 - Pos22*Pos13;
    const double Adj23 = Pos12*Pos13 - Pos11*Pos23;
    const double Adj33 = Pos11*Pos22 - Pos12*Pos12;

#ifdef USE_INVARIANT_POS_IND
    // extra computation so as to use matrix invariants
    const double Adj11 = Pos22*Pos33 - Pos23*Pos23;
    const double Adj22 = Pos33*Pos11 - Pos13*Pos13;

#endif /* USE_INVARIANT_POS_IND */
    const double detPos = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;

    // positive-definite means that all three of these are positive.
    //if(Pos11<0.||Adj33<0.||detPos<0.) retval = false;
#ifndef USE_INVARIANT_POS_IND
    pos.set(mp,1, Pos11);
    pos.set(mp,2, Adj33);
#else /* USE_INVARIANT_POS_IND */
    pos.set(mp,1, Pos11+Pos22+Pos33);
    pos.set(mp,2, Adj11+Adj22+Adj33);
#endif /* USE_INVARIANT_POS_IND */
    pos.set(mp,3, detPos);
    // is this the worst offender so far?
    for(int i=1;i<=3;i++)
    {
      // ideally the worst offender would simply 
      // be the positivity tensor with the lowest minimal eigenvalue,
      // but that is expensive to compute.
      const double value = pos.get(mp,i);
      const double offense_severity = -value/pos0[i];
      if(offense_severity > max_offense_severity)
      {
        max_offense_severity = offense_severity;
        retval.set(mp,i,offense_severity,value);
      }
    }
  }
  if(retval.is_set())
  {
    if(debug3){ printf("worst_offender: "); retval.print(); printf("\n"); }
  }
  assert_gt(retval.get_poly_idx(),0);
  return retval;
}

// smallest_positive_zero(poly3) <=
// smallest_positive_zero(poly2) <=
// smallest_positive_zero(poly1),
// so we really only need to get poly3
// if we aim for the smallest positive zero.
// But in practice we are satisfied with any
// zero that makes all positivity indicators positive.
// (issue: unless we choose the smallest eigenvalue
// the recipe is not smooth.)
void get_gas10_positivity_polynomials_at_point(
  int mp,
  Polynomial& poly1,
  Polynomial& poly2,
  Polynomial& poly3,
  const dTensor1& Q0,
  const dTensor2& dQ,
  double epsilon)
{
  using namespace TenMomentComponentID;
  // conserved variables are linear polynomials.
  //
  const Polynomial rho(Q0.get(_rho),dQ.get(mp,_rho));
  const Polynomial M1(Q0.get(_M1),dQ.get(mp,_M1));
  const Polynomial M2(Q0.get(_M2),dQ.get(mp,_M2));
  const Polynomial M3(Q0.get(_M3),dQ.get(mp,_M3));
  //
  const Polynomial Ne11(Q0.get(_N11) - epsilon, dQ.get(mp,_N11));
  const Polynomial Ne12(Q0.get(_N12)          , dQ.get(mp,_N12));
  const Polynomial Ne13(Q0.get(_N13)          , dQ.get(mp,_N13));
  const Polynomial Ne22(Q0.get(_N22) - epsilon, dQ.get(mp,_N22));
  const Polynomial Ne23(Q0.get(_N23)          , dQ.get(mp,_N23));
  const Polynomial Ne33(Q0.get(_N33) - epsilon, dQ.get(mp,_N33));

  // quadratic polynomials
  // 
#ifndef USE_INVARIANT_POS_IND
  Polynomial& Pos11 = poly1;
                   Pos11 = rho*Ne11 - M1*M1;
#else /* USE_INVARIANT_POS_IND */
  const Polynomial Pos11 = rho*Ne11 - M1*M1;
#endif /* USE_INVARIANT_POS_IND */
  const Polynomial Pos12 = rho*Ne12 - M1*M2;
  const Polynomial Pos13 = rho*Ne13 - M1*M3;
  const Polynomial Pos22 = rho*Ne22 - M2*M2;
  const Polynomial Pos23 = rho*Ne23 - M2*M3;
  const Polynomial Pos33 = rho*Ne33 - M3*M3;

  // quartic polynomials
#ifndef USE_INVARIANT_POS_IND
  Polynomial& Adj33 = poly2;
                   Adj33 = Pos11*Pos22 - Pos12*Pos12;
#else /* USE_INVARIANT_POS_IND */
  const Polynomial Adj33 = Pos11*Pos22 - Pos12*Pos12;
#endif /* USE_INVARIANT_POS_IND */
  const Polynomial Adj13 = Pos12*Pos23 - Pos22*Pos13;
  const Polynomial Adj23 = Pos12*Pos13 - Pos11*Pos23;
#ifdef USE_INVARIANT_POS_IND
  // extra computation so as to use matrix invariants
  const Polynomial Adj11 = Pos22*Pos33 - Pos23*Pos23;
  const Polynomial Adj22 = Pos33*Pos11 - Pos13*Pos13;
#endif /* USE_INVARIANT_POS_IND */
  // det(P) (hexic polynomial)
#ifdef USE_INVARIANT_POS_IND
  poly1 = Pos11+Pos22+Pos33;
  poly2 = Adj11 + Adj22 + Adj33;
#endif /* USE_INVARIANT_POS_IND */
  poly3 = Adj13*Pos13 + Adj23*Pos23 + Adj33*Pos33;

  // the root-finder assumes fixed machine-epsilon
#ifndef USE_INVARIANT_POS_IND
  poly2; //*= 1e16;
  poly3; //*= 1e32;
#else /* USE_INVARIANT_POS_IND */
  //poly2; //*= 1e16;
  //poly3; //*= 1e32;
#endif /* USE_INVARIANT_POS_IND */
}

class PolynomialArray2
{
  Polynomial*vec;
  int s1;
  int s2;
  int s;
 private:
  void operator=(const PolynomialArray2&); // disable assignment
  PolynomialArray2(const PolynomialArray2&); // disable copy constructor
  const int get_k(int n1, int n2) const
  {
      int k = (n1-1)*s2 + (n2-1);
      assert_gt(n1,0); assert_le(n1,s1);
      assert_gt(n2,0); assert_le(n2,s2);
      return k;
  }
 public:
  PolynomialArray2(int _s1, int _s2):s1(_s1),s2(_s2)
  {
    s = s1*s2;
    vec = new Polynomial[s];
  }
  ~PolynomialArray2()
  {
    delete [] vec;
  }
  Polynomial& ref(int n1, int n2)
  {
      return vec[get_k(n1,n2)];
  }
  void set(int n1, int n2, const Polynomial& p)
  {
      vec[get_k(n1,n2)] = p;
  }
};

// enforce positivity at all points
//
// Let q(t) = q0 + t*(q-q0), where q0 is cell average.
//
// At each (positivity) point m:
//   Pos11:=pos_1:=pos(m,1) is a quadratic polynomial in t,
//   Adj33:=pos_2:=pos(m,2) is a fourth-order polynomial in t,
//  detPos:=pos_3:=pos(m,3) is a sixth-order polynomial in t.
//
// Find a theta for which one of the positivity indicators 
// of q(theta) is zero (probably detPos) and the others are positive.
//
double Gas10_CellPositivityLimiter::get_theta_using_pos_indicators(
  Offender& worst_offender, dTensor2& pos,
  const dTensor1& Q0,const dTensor2& dQ,
  double epsilon)const
{
  //const int numPoints = Legendre2d::get_numPositivityPoints();
  const int numPoints = points.get_numPoints();
  const int numPolynomials = 3*numPoints;
  PolynomialArray2 parr(numPoints,3);
  iTensor1 pset(numPoints);
  pset.setall(0);

  // Try a shortcut: find the worst-offending polynomial
  // and find one of its roots between 0 and 1
  // and see if that value makes everything else positive.
  //
  // hopefully this loop will almost always be executed only once
  // and very rarely many times.
  // (if executed many times maybe it would be more efficient to compute
  // all the polynomials and work with sturm sequences.)
  //
  double theta_upper = 1.;
  for(int ctr=0; ctr<= 12*numPoints; ctr++)
  {
    const int mp = worst_offender.get_point_idx();
    if(!pset.get(mp))
    {
      // Polynomial poly[4];
      get_gas10_positivity_polynomials_at_point(mp,
        parr.ref(mp,1),
        parr.ref(mp,2),
        parr.ref(mp,3),
        Q0, dQ, epsilon);
      // otherwise numerical error can cause code to crash
      // (if doing this then return 0 if value is negative
      // even for cell average.)
      //if(LimiterParams::get_pos1_shift())
      //  parr.ref(mp,1).subtract(LimiterParams::get_pos1_shift());
      //if(LimiterParams::get_pos2_shift())
      //  parr.ref(mp,2).subtract(LimiterParams::get_pos2_shift());
      //if(LimiterParams::get_pos3_shift())
      //  parr.ref(mp,3).subtract(LimiterParams::get_pos3_shift());
      pset.set(mp,1);
    }
    // Find a positive root of one of the indicator polynomials that 
    // makes the other polynomials positive.
    // (note that dividing by (root-x) does not change whether
    // the polynomial is positive for values less than root.)
    double try_theta = 0.;
    int pidx = worst_offender.get_poly_idx();
    dprintf3("\n  parr.ref(%d,%d).eval(%24.16e))=%24.16e",
      mp,pidx,theta_upper, parr.ref(mp,pidx).eval(theta_upper));
    parr.ref(mp,pidx).isolate_and_divide_out_root_minus_x(try_theta,theta_upper);
    dprintf3("\n  parr.ref(%d,%d).eval(%24.16e))=%24.16e",
      mp,pidx,try_theta, parr.ref(mp,pidx).eval(try_theta));
    //assert_ge(parr.ref(mp,pidx).eval(try_theta),0.);
    assert_le(try_theta,theta_upper);

    // cycle through polynomials until all are positive at try_theta;
    // there are potentially 2+4+6=12 roots associated with this point.
    int k;
    for(k=1;k<=13;k++)
    {
      for(int j=3;j>=1;j--)
      {
        if(j!=pidx) // (the indicator for this one is currently zero)
        {
          const double value = parr.ref(mp,j).eval(try_theta);
          if(value < 0.)
          {
            // every time we change theta we test the other indicators again
            pidx = j;
            theta_upper = try_theta;
            try_theta = 0.;
            dprintf2("\n  parr.ref(%d,%d).eval(%24.16e))=%24.16e",
              mp,j,theta_upper, parr.ref(mp,j).eval(theta_upper));
            parr.ref(mp,j).isolate_and_divide_out_root_minus_x(try_theta,theta_upper);
            dprintf3("\n  parr.ref(%d,%d).eval(%24.16e))=%24.16e",
              mp,j,try_theta, parr.ref(mp,j).eval(try_theta));
            goto try_again;
          }
        }
        goto all_3_indicators_are_nonnegative;
      }
      try_again: ;
    }
    all_3_indicators_are_nonnegative:
    assert_le(k,12);
    //double try_theta = poly3.get_smallest_positive_root();

    // are the (other) points now okay?
    Offender old_worst_offender = worst_offender;
    worst_offender = get_gas10_positivity_indicators_at_points(
      pos,Q0,dQ,try_theta,epsilon);
    if(worst_offender.not_set())
    {
      dprintf3("i=%d,j=%d,moffset=%d,k=%d",i,j,moffset,k);
      return try_theta;
    }
    else
    {
      theta_upper = try_theta;
    }
    // repeated worst offender?
    if(old_worst_offender.get_point_idx()==worst_offender.get_point_idx())
    {
      // presumably we are here because of numerical error.
      // by assumption the cell average is okay.
      // so it *is* possible to decrease theta sufficiently to 
      // satisfy positivity.  Shall we try to do that somehow?

      if(debug2) // why is the worst offender still set?
      {
        println("worst offender still set.")
        printvn(epsilon);
        print_evals_and_pos_indicators_for_point(mp,Q0,dQ,pos);
      }
      return try_theta;

      // save offender data
      const int poly_idx = worst_offender.get_poly_idx();
      const double pvt = pos.get(worst_offender.get_point_idx(),poly_idx);
      dprint(pvt);
      if(debug2){
        printf("worst_offender: ");
        worst_offender.print(); printf("\n");
      }

      // do we satisfy a softer requirement?
      //
      const double soft_epsilon = epsilon*.75;
      worst_offender = get_gas10_positivity_indicators_at_points(
        pos,Q0,dQ,try_theta,soft_epsilon);
      if(debug2)
      {
        printf("i=%d,j=%d,moffset=%d,k=%d (soft-epsilon),"
               "\n  try_theta=%24.16e, worst_offender: ",
               i,j,moffset,k,try_theta);
        worst_offender.print();
        printf("\n");
      }
      if(worst_offender.not_set())
      {
        return try_theta;
      }

      // use offender data to attempt linear interpolation
      if(false) // this does not help
      {
        double pos00[4];
        assert(get_gas10_positivity_indicators_for_cell_avg_state(
          Q0,pos00,epsilon));
        const double pv0 = pos00[poly_idx];
        dprint(pv0);
        const double lin_theta = try_theta*(pv0/(pv0-pvt));

        // does the linear interpolation satisfy the softer requirement?
        //
        if(lin_theta>0. && lin_theta < try_theta)
        {
          worst_offender = get_gas10_positivity_indicators_at_points(
            pos,Q0,dQ,lin_theta,soft_epsilon);
          if(debug2)
          {
            printf("i=%d,j=%d,moffset=%d,k=%d (soft-epsilon, linear),"
                   "\n  lin_theta=%24.16e, worst_offender=",
                   i,j,moffset,k,lin_theta);
            worst_offender.print(); printf("\n");
          }
          if(worst_offender.not_set())
          {
            return lin_theta;
          }
        }
      }

      // give up and set theta=0
      //
      if(false)
      {
        dprintf("i=%d,j=%d,moffset=%d,k=%d"
                "\n  zeroing repeated offender: ",
                i,j,moffset,k);
        worst_offender.print(); printf("\n");
      }
    }
  }

  eprintf("failed to find theta that satisfies positivity everywhere\n");
}

// this function is a cubic spline which satisfies
// f(0) = max_v,
// f(max_x) = 0,
// f'(0) = f'(max_x) = 0.
inline double cubic_transition(double x, double max_x, double max_v)
{
  if(x<=0.) return max_v;
  if(x>=max_x) return 0.;
  const double x_m = x/max_x;
  return max_v*(1 - x_m)*(1 - x_m)*(1+2*x_m);
  //return max_v*(1./(max_x*max_x*max_x))*(max_x - x)*(max_x - x)*(max_x+2*x);
}

// return true if not modifying solution
bool gas10_enforce_positivity_of_cell_average(
    int i,int j,int moffset, dTensorBC4& q, double* pos0)
{
  // check positivity of cell average
  dTensor1 Q0(10);
  ::get_Q0_in_cell(i,j,moffset,q,Q0);
  const double avg_eps = LimiterParams::get_avg_pressure_epsilon(); // 0.;
  const bool ispos = get_gas10_positivity_indicators_for_cell_avg_state(
    Q0, pos0, avg_eps);

  double iso_factor = 1.; // initialization
  // isotropize at most this amount unless positivity is violated
  // (and enforce that pressure anisotropy is no greater
  // than about a factor of 1/max_iso)
  const double max_iso = 1e-4;
  if(!ispos)
  {
    double evals[3];
    dTensor1 prim(10);
    gas10_convert_cons_to_prim(0, Q0, prim);
    double p = gas10_get_pressure_evals_from_prim(prim, evals);
    if(debug2)
    {
      printf(" for cell:");
      printvc(i);
      printvc(j);
      printvc(moffset);
      println("\n  isotropizing pressure evals:")
      printvn(evals[0]);
      printvn(evals[1]);
      printvn(evals[2]);
      println("  because of positivity indicators:")
      printvn(pos0[0]);
      printvn(pos0[1]);
      printvn(pos0[2]);
      printvn(pos0[3]);
    }
    if(evals[0]>avg_eps)
    {
      if(debug2){
        printf("positivity of cell avg seems not actually a problem:");
        printv(evals[0]);
      }
      iso_factor = (1.-max_iso);
    }
    else
    {
      // isotropize just enough so that (ignoring max_iso) the
      // minimal eigenvalue becomes avg_eps
      iso_factor = (1.-max_iso)*(p-avg_eps)/(p-evals[0]);
    }
  }
  else // for smooth transition to isotropization
  {
    assert(ispos);
    // estimate ratio of smallest to largest eigenvalue
    //             cases:               ellips prolate oblate   isotropic
    const double& trPos = pos0[1]; // ~ L3     L3      2*L3     3*L3
    const double& trAdj = pos0[2]; // ~ L2*L3  2*L2*L3 L2*L3    3*L2*L3
    const double& detPos= pos0[3]; // = L1*L2*L3 (exactly)
    // for smoothness isotropization we want to base isotropization factor
    // on the ratio L1_L3.
    // this estimate is between L1/L3 and (1/9)*(L1/L3);
    // for small L1_L3 it is between L1/L3 and about (1/2)*(L1/L3):
    // ~(1/9)*(L1/L3) ~= 1/9 in the isotropic case L1 ~= L2 ~= L3,
    // ~(1/2)*(L1/L3) in the severely prolate case L1 ~= L2 << L3,
    // ~(1/2)*(L1/L3) in the severely oblate  case L1 << L2 ~= L3,
    // ~      (L1/L3) in the very ellipsoidal case L1 << L2 << L3.
    const double L1_L3_estimate = detPos/(trPos*trAdj);
    // could instead isotropize to enforce that this is greater
    // than epsilon:
    //const double L1_estimate = detPos/trAdj;
    //
    // smoothly transition to corrective isotropization
    //
    const double trigger = 2.5*max_iso;
    if(L1_L3_estimate < trigger)
    {
      assert_gt(L1_L3_estimate,0.);
      const double iso = cubic_transition(L1_L3_estimate, trigger, max_iso);
      iso_factor = 1. - iso;
      if(debug3)
      {
        double evals[3];
        dTensor1 prim(10);
        gas10_convert_cons_to_prim(0, Q0, prim);
        double p = gas10_get_pressure_evals_from_prim(prim, evals);
        printf(" for cell:");
        printvc(i);
        printvc(j);
        printvc(moffset);
        println("\n  isotropizing pressure evals:")
        printvn(evals[0]);
        printvn(evals[1]);
        printvn(evals[2]);
        println("  because of pressure anisotropy:")
        printvn(L1_L3_estimate);
        println("  positivity indicators are:")
        printvn(pos0[0]);
        printvn(pos0[1]);
        printvn(pos0[2]);
        printvn(pos0[3]);
      }
    }
  }

  assert_le(iso_factor,1.);
  assert_ge(iso_factor,0.);
  const bool modify_solution = (iso_factor < 1.);
  if(modify_solution) // modify solution
  {
    dTensor1 prim(10);
    gas10_convert_cons_to_prim(0, Q0, prim);
    gas10_rescale_pressure_anisotropy(prim, iso_factor);
    gas10_convert_prim_to_cons(0, prim, Q0);

    const bool check_pos = true;
    if(check_pos)
    {
      double iso_pos[4];
      const bool iso_ispos = get_gas10_positivity_indicators_for_cell_avg_state(
        Q0, iso_pos, 0.); //LimiterParams::get_verify_avg_pressure_epsilon());
      if(!iso_ispos)
      {
        eprintf("cell average still violates positivity of pressure");
      }
      if(debug3)
      {
        double iso_evals[3];
        gas10_get_pressure_evals_from_prim(prim, iso_evals);
        println(" isotropized pressure evals:");
        printvn(iso_evals[0]);
        printvn(iso_evals[1]);
        printvn(iso_evals[2]);
        println("  and new positivity indicators:")
        printvn(iso_pos[0]);
        printvn(iso_pos[1]);
        printvn(iso_pos[2]);
        printvn(iso_pos[3]);
      }
    }
    ::put_Q0_in_cell(i,j,moffset,q,Q0);
  }
  return !modify_solution;
}

bool gas10_verify_positivity_of_cell_average(
    int i,int j,int moffset, const dTensorBC4& q, double* pos0)
{
  // check positivity of cell average
  dTensor1 Q0(10);
  ::get_Q0_in_cell(i,j,moffset,q,Q0);
  return get_gas10_positivity_indicators_for_cell_avg_state(Q0, pos0, 0.);
}

// return true if enforcement was unnecessary
bool gas10_enforce_positivity_of_cell_averages(
  int moffset, dTensorBC4& q, double* pos0min)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mbc  = q.getmbc();

  bool all_okay = true;
  // loop through cells
  //
  // parallelization needs intermediate "accumulators"
  dTensor2d minPos0(mx,4,1,0);
  #pragma omp parallel for
  for (int i=1; i<=mx; i++)   
  {
    // initialize inner "accumulator" for minimum positivity indicators
    double pos0min[4]; for(int k=0;k<4;k++) pos0min[k]=DBL_MAX;
    for (int j=1; j<=my; j++)   
    {
      double pos0[4];
      bool okay = gas10_enforce_positivity_of_cell_average(i,j,moffset,q, pos0);
      all_okay = all_okay && okay;
      if(!okay&&debug3){
        printvc(i);
        printvc(j);
        printf("\n");
        printvn(pos0[0])
        printvn(pos0[1])
        printvn(pos0[2])
        printvn(pos0[3])
      }
      for(int k=0;k<4;k++) pos0min[k]=Min(pos0[k],pos0min[k]);
    }
    // save minimum positivity values for this loop
    for(int k=0;k<4;k++) minPos0.set(i,k, pos0min[k]);
  }
  // initialize outer "accumulator" for minimum positivity indicators
  for(int k=0;k<4;k++) pos0min[k]=DBL_MAX;
  for (int i=1; i<=mx; i++)   
  {
    for(int k=0;k<4;k++) pos0min[k]=Min(minPos0.get(i,k),pos0min[k]);
  }
  //for(int k=0;k<4;k++) if(pos0min[k]<0) return false;
  //return true;
  return all_okay;
}

bool gas10_verify_positivity_of_cell_averages(
  int moffset, const dTensorBC4& q, double* pos0min)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mbc  = q.getmbc();

  // loop through cells
  //
  // parallelization needs intermediate "accumulators"
  dTensor2d minPos0(mx,4,1,0);
  #pragma omp parallel for
  for (int i=1; i<=mx; i++)   
  {
    // initialize inner "accumulator" for minimum positivity indicators
    double pos0min[4]; for(int k=0;k<4;k++) pos0min[k]=DBL_MAX;
    for (int j=1; j<=my; j++)   
    {
      //ApplyPositivityLimiterInCell(i,j,q);
      double pos0[4];
      bool ispos = gas10_verify_positivity_of_cell_average(i,j,moffset,q, pos0);
      if(!ispos && debug1)
      // print evals for cell avg
      {
        double evals[3];
        dTensor1 Q0(10);
        get_Q0_in_cell(i,j,moffset,q, Q0);
        gas10_get_pressure_evals(Q0, evals);
        printf(" for cell:");
        printvc(i);
        printvc(j);
        printvc(moffset);
        println("\n  pressure evals:")
        printvn(evals[0]);
        printvn(evals[1]);
        printvn(evals[2]);
        println("  positivity indicators:")
        printvn(pos0[0]);
        printvn(pos0[1]);
        printvn(pos0[2]);
        printvn(pos0[3]);
      }
      for(int k=0;k<4;k++) pos0min[k]=Min(pos0[k],pos0min[k]);
    }
    // save minimum positivity values for this loop
    for(int k=0;k<4;k++) minPos0.set(i,k, pos0min[k]);
  }
  // initialize outer "accumulator" for minimum positivity indicators
  for(int k=0;k<4;k++) pos0min[k]=DBL_MAX;
  for (int i=1; i<=mx; i++)   
  {
    for(int k=0;k<4;k++) pos0min[k]=Min(minPos0.get(i,k),pos0min[k]);
  }

  // rescale pressure positivity indicators to make them comparable?
  // I'm not totally clear what the rescaling should be. The
  // positivity indicators are smooth (polynomial) functions of
  // the state. The commented-out non-smooth rescaling below
  // cannot be right. In case all three eigenvalues go to
  // zero simultaneously P11 would be linear and would be the
  // dominating negative eigenvalue. In case only one eigenvalue
  // goes to zero we can expect the indicators to have nonzero
  // derivatives as they cross zero.
  //
  //double p2 = pos0min[2];
  //double p3 = pos0min[3];
  //pos0min[2] = std::signbit(pos0min[2])*sqrt(fabs(pos0min[2]));
  //pos0min[3] = std::signbit(pos0min[3])*cbrt(fabs(pos0min[3]));

  for(int k=0;k<4;k++) if(pos0min[k]<0) return false;
  return true;
}

int Gas10_CellPositivityLimiter::check()const
{
  const double pressure_epsilon = LimiterParams::get_pressure_epsilon();

  using namespace PositivityReturnType;
  int ret = MAX;

  ret = iMin(ret, checkDensity());
  // not much point in checking the pressure if the density is negative
  if(ret == NEG_DENSITY) return ret;

  IntervalArray Qinterval(11);
  get_extreme_possible_states(10,Qinterval);
  if(gas10_guaranteed_positivity_using_intervals(
       Qinterval,pressure_epsilon))
  {
    dprintf3("for (%d,%d) INTERVALS_CONFIRMED_POSITIVITY",i,j);
    ret = iMin(ret,INTERVALS_CONFIRMED_POSITIVITY);
    return ret;
  }

  // check positivity of cell average
  dTensor1 Q0(10);
  get_Q0_in_cell(Q0);
  double pos0[4];
  bool okay = get_gas10_positivity_indicators_for_cell_avg_state(
    Q0,pos0,pressure_epsilon);
  assert_gt(pos0[0], 0.); // already checked this
  if(!okay)
  {
    dprint1(pos0[0]);
    dprint1(pos0[1]);
    dprint1(pos0[2]);
    dprint1(pos0[3]);
    dprintf1("cell average violates positivity of pressure");
    ret = iMin(ret, NEG_PRESSURE);
    return ret;
  }

  const int numPoints = points.get_numPoints();
  dTensor1& Q = Q0;
  for(int mp=1;mp<=numPoints;mp++)
  {
    get_Q_at_point(Q,mp);

    int gas10_check_positivity_cons(int moffset, const dTensor1& Q);
    int ret_mp = gas10_check_positivity_cons(0, Q);
    ret = iMin(ret,ret_mp);
    if(ret_mp<0)
    {
      dprintf("positivity violation for:"
        "\n  moffset = %d,"
        "\n  mp = %d,"
        "\n  i = %d,"
        "\n  j = %d,",
        moffset, mp, i, j);
    }
  }
  return ret;
}

void Gas10_CellPositivityLimiter::print_evals_and_pos_indicators_for_point(
  int mp,const dTensor1& Q0,const dTensor2& dQ,const dTensor2& pos)const
{
  dTensor1 Q(10);
  ::get_Q_at_point(Q, mp,Q0,dQ);
  double evals[3];
  gas10_get_pressure_evals(Q, evals);
  printf(" at point:");
  printvc(moffset);
  printvc(mp);
  printvc(i);
  printv(j);
  println("\n  pressure evals:")
  printvn(evals[0]);
  printvn(evals[1]);
  printvn(evals[2]);
  println("  positivity indicators:")
  printvn(pos.get(mp,1));
  printvn(pos.get(mp,2));
  printvn(pos.get(mp,3));
  double eval0[3];
  gas10_get_pressure_evals(Q0, eval0);
  println(" pressure evals at cell average state:")
  printvn(eval0[0]);
  printvn(eval0[1]);
  printvn(eval0[2]);
  double pos00[4];
  assert(get_gas10_positivity_indicators_for_cell_avg_state(
    Q0,pos00,LimiterParams::get_verify_avg_pressure_epsilon()));
  println(" positivity indicators for cell average state:");
  printvn(pos00[0]);
  printvn(pos00[1]);
  printvn(pos00[2]);
  printvn(pos00[3]);
}

int Gas10_CellPositivityLimiter::apply()
{
  int ret = PositivityReturnType::MAX;
  // enforce positivity of density
  double min_density;
  ret = iMin(ret,applyDensityLimiter(&min_density));

  // remainder enforces positivity of pressure tensor

  const double pressure_epsilon = LimiterParams::get_pressure_epsilon();
  const double avg_eps = LimiterParams::get_avg_pressure_epsilon();

  // can we confirm positivity merely using magnitudes
  // of Legendre coefficients and interval arithmetic?
  // (for refined mesh the lower-order terms (in particular
  // the terms of the linear part) should dominate,
  // so this should usually be enough to confirm positivity
  //
  IntervalArray Qinterval(11);
  get_extreme_possible_states(10,Qinterval);
  // since positivity-limiting has been applied to the density we in fact have:
  Qinterval.fetch(1).set_min(min_density);
  if(gas10_guaranteed_positivity_using_intervals(Qinterval,avg_eps))
  {
    dprintf3("for (%d,%d) INTERVALS_CONFIRMED_POSITIVITY",i,j);
    ret = iMin(ret, INTERVALS_CONFIRMED_POSITIVITY);
    return ret;
  }

  // check positivity of cell average
  dTensor1 Q0(10);
  get_Q0_in_cell(Q0);
  double pos0[4];
  bool okay = get_gas10_positivity_indicators_for_cell_avg_state(
    Q0,pos0,avg_eps);
  if(!okay)
  {
    double pos00[4];
    bool okay0 = get_gas10_positivity_indicators_for_cell_avg_state(
      Q0,pos00,0.);
    if(!okay0)
    {
      dprint1(pos00[0]);
      dprint1(pos00[1]);
      dprint1(pos00[2]);
      dprint1(pos00[3]);
      eprintf("positivity violation:");
      //timeStepException.Throw(NEGATIVE_PRESSURE);
    }
    // Compute the minimal eigenvalue of the pressure tensor and
    // shift the pressure tensor so that the minimal eigenvalue
    // of P equals avg_eps. This is expensive; we really just
    // need to do something smooth that makes the pressure exceed
    // epsilon by a certain amount but not more than a certain
    // greater amount (based, e.g., on the positivity indicator
    // values).
    double evals[3];
    gas10_get_pressure_evals(Q0, evals);
    const double pmin = evals[0];
    assert_ge(pmin,0.);
    assert_le(pmin,avg_eps);
    const double pshift = avg_eps-pmin;
    dprintf("shifting pressure by %24.16e",pshift);
    shift_pressure(Q0,pshift);
  }
  // do_we_have_positivity_at_all_the_positivity_points?
  // (intervals were not enough, so check state at positivity points.)
  //
  const int numPoints = points.get_numPoints();
  //const int numPoints = Legendre2d::get_numPositivityPoints();
  dTensor2 dQ(numPoints,10);
  get_dQ_at_positivity_points(dQ);
  // check pressure tensor positivity at all points
  dTensor2 pos(numPoints,3);
  Offender worst_offender
    = get_gas10_positivity_indicators_at_points(pos,Q0,dQ,1.,pressure_epsilon);
  if(worst_offender.not_set())
  {
    //dprintf("for (%d,%d,%d) EVALUATIONS_CONFIRMED_POSITIVITY",i,j,moffset);
    ret = iMin(ret, EVALUATIONS_CONFIRMED_POSITIVITY);
    return ret;
  }

  dprintf3("for (%d,%d,%d) worst_offender: %d.%d: %24.16e",
    i,j,moffset,
    worst_offender.get_point_idx(),
    worst_offender.get_poly_idx(),
    worst_offender.get_severity());
  
  const double theta
    = get_theta_using_pos_indicators(worst_offender,pos,Q0,dQ,pressure_epsilon);
  rescale_dq_in_cell(i,j,moffset,10,q,theta); // assign q := q0+theta(q-q0)
  // print out minimal pressure eigenvalue

  // confirm positivity at all points
  // (take out this expensive and redundant check when debugged)
  //
  if(true)
  {
    get_dQ_at_positivity_points(dQ); // update dQ from q
    worst_offender = get_gas10_positivity_indicators_at_points(
        pos,Q0,dQ,1.,LimiterParams::get_verify_pressure_epsilon());
    if(worst_offender.is_set())
    {
      if(debug1) {
        printf("failed to enforce positivity using theta = %24.16e;\n"
          "  worst_offender = ",theta);
        worst_offender.print(); printf("\n");
      }
    }
    // show eigenvalues for worst-offending point
    if(debug3)
    {
      int mp = worst_offender.get_point_idx();
      print_evals_and_pos_indicators_for_point(mp,Q0,dQ,pos);
    }
  }
  dprintf3("for (%d,%d,%d) ENFORCED_POSITIVITY",i,j,moffset);
  ret = iMin(ret, ENFORCED_PRESSURE);
  return ret;
}

// =============== full model stuff ================

// check positivity in cell
// bail as soon as we can.
//
// as the mesh is refined higher-order coefficients should
// become small.  So for most cells we can quickly confirm
// positivity using the magnitudes of higher coefficients.
//
int PositivityLimiter::checkPositivityInCell(int i,int j) const
{
  using namespace PlasmaModelType;
  int ret = PositivityReturnType::MAX;
  int Gas10offset = 10;
  if(plasmaParams.get_model()==g20)
    Gas10offset = 20;
  switch(plasmaParams.get_model())
  {
   case g20:
   case g10:
   {
    Gas10_CellPositivityLimiter elc_limiter(i,j,Gas10offset,q,points);
    ret = iMin(ret,elc_limiter.check());
   }
   case p20:
   case p10:
   {
    Gas10_CellPositivityLimiter ion_limiter(i,j,0,q,points);
    ret = iMin(ret,ion_limiter.check());
   }
    break;
   case g05:
   {
    Gas05_CellPositivityLimiter elc_limiter(i,j,5,q,points);
    ret=iMin(ret,elc_limiter.check());
   }
   case p05:
   case mhd: // good enough?
   {
    Gas05_CellPositivityLimiter ion_limiter(i,j,5,q,points);
    ret = iMin(ret,ion_limiter.check());
   }
    break;
   case i10e5:
   {
    Gas10_CellPositivityLimiter ion_limiter(i,j,0,q,points);
    ret = iMin(ret,ion_limiter.check());
    Gas05_CellPositivityLimiter elc_limiter(i,j,10,q,points);
    ret = iMin(ret,elc_limiter.check());
   }
    break;
   default:
    unsupported_value_error(plasmaParams.get_model());
  }
  return (PositivityReturnType::Enum) ret;
}

// check if positivity limiters are needed.
// bail as soon as we can.
//
// as the mesh is refined higher-order coefficients should
// become small.  So for most cells we can quickly confirm
// positivity using the magnitudes of higher coefficients.
//
int PositivityLimiter::ApplyInCell(int i, int j)
{
  // compute lower and upper bounds on each component
  // (based on values of Legendre polynomials within the cell)
  using namespace PlasmaModelType;
  int ret = PositivityReturnType::MAX;
  int Gas10offset = 10;
  if(plasmaParams.get_model()==g20)
    Gas10offset = 20;
  // assert(plasmaParams.get_model()==PlasmaModelType::i10e5);
  switch(plasmaParams.get_model())
  {
   case g20:
   case g10:
   {
    Gas10_CellPositivityLimiter elc_limiter(i,j,Gas10offset,q,points);
    ret = iMin(ret,elc_limiter.apply());
   }
   case p20:
   case p10:
   {
    Gas10_CellPositivityLimiter ion_limiter(i,j,0,q,points);
    ret = iMin(ret,ion_limiter.apply());
   }
    break;
   case g05:
   {
    Gas05_CellPositivityLimiter elc_limiter(i,j,5,q,points);
    ret=iMin(ret,elc_limiter.apply());
   }
   case p05:
   case mhd: // good enough?
   {
    Gas05_CellPositivityLimiter ion_limiter(i,j,0,q,points);
    ret = iMin(ret,ion_limiter.apply());
   }
    break;
   case i10e5:
   {
    Gas10_CellPositivityLimiter ion_limiter(i,j,0,q,points);
    ret = iMin(ret,ion_limiter.apply());
    Gas05_CellPositivityLimiter elc_limiter(i,j,10,q,points);
    ret = iMin(ret,elc_limiter.apply());
   }
    break;
   default:
    unsupported_value_error(plasmaParams.get_model());
  }
  return ret;
}

#include "Legendre2d.h"

// verify positivity at positivity points
//
int PositivityLimiter::checkPositivity()const
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mbc  = q.getmbc();

  // loop through cells
  //
  // need to check positivity in neighbors of cells that will be evolved
  int ret = PositivityReturnType::MAX;
  #pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)   
  {
    for (int j=(1-mbc); j<=(my+mbc); j++)   
    {
      ret = iMin(ret,checkPositivityInCell(i,j));
    }
  }
  return ret;
}

// apply positivity limiter to guarantee physicality at a set a points
//
// + positivity at guass-lobatto points provably is maintained
//   for sufficiently short time step. Eigenvalues (speeds) are
//   computed at guassian quadrature points on the boundaries
//   (which by design coincide with gauss-lobatto points). For
//   limiting, eigenstructure is computed for the cell average.
//
int PositivityLimiter::apply()
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mbc  = q.getmbc();

  // loop through cells
  //
  // need to have positivity in neighbors of cells that will be evolved;
  // might as well have it in all.
  //
  // any exception thrown within a thread must be caught in the thread
  int ret = PositivityReturnType::MAX;
  #pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)   
  {
    for (int j=(1-mbc); j<=(my+mbc); j++)   
    {
      ret = iMin(ret,ApplyInCell(i,j));
    }
  }
  if(ret <= PositivityReturnType::ENFORCED_POSITIVITY)
  {
    dprintf("ret: %s",get_string((PositivityReturnType::Enum) ret));
  }
  return ret;
}
PositivityPointInfo::PositivityPointInfo(
  PositivityPointType::Enum pointType_in):
  pointType(pointType_in),
  numPoints          (0),
  polynomValsAtPoints(0),
  polynomIntervals   (0)
{
  set_polynomIntervals(&Legendre2d::get_phi_interval());
  using namespace PositivityPointType;
  switch(pointType)
  {
    case FLUX:
      set_polynomValsAtPoints(&Legendre2d::get_phiAtPositivityPoints());
      break;
     
    case SOURCE:
      set_polynomValsAtPoints(&Legendre2d::get_phi());
      break;

    case RIEMANN:
      set_polynomValsAtPoints(&Legendre2d::get_phiAtRiemannPoints());
      break;

    default:
      unsupported_value_error(pointType);
  }
  numPoints = polynomValsAtPoints->getsize(1);
}

// *** end positivity_limiters ***
