#ifndef Relax_h
#define Relax_h
#include <cmath>
#include "assert.h"
#include "gas10.h"
#include "tensors.h"
#include "PlasmaParams.h"
#include "float.h"

// in-line methods for relaxation
//
// (inline makes sense for short functions
// with many arguments that are called infrequently
// in the code and frequently during execution)
//
namespace Relax{

inline double get_portion_to_keep(double iso_rate,double d_time)
{
    double portion_to_keep;
    if(iso_rate <= DBL_MAX)
    {
        const double decay = d_time*iso_rate;
        // RKn approximates exp(-decay) with its n-th order
        // Taylor series approximation "exp".  Setting
        // "exp"(-decay)=1 and using the formula for
        // the roots of a cubic reveals that for RK3
        // max_decay = 1/r+1-r =~= 2.5127453266183282,
        // where r=cbrt(sqrt(17)-4) =~= 0.49746129930366850
        dprintf3("decay cfl = %24.16e",decay*(1./2.5127453266183282));
        portion_to_keep = exp(-decay);
    }
    else
        portion_to_keep = 0.;
    return portion_to_keep;
}

inline double get_iso_rate(
  IsoPeriodType::Enum iso_period_type, double in_iso_rate, double pcl_mass_inv,
  double rho, double p,
  double P11, double P12, double P13, double P22, double P23, double P33)
{
  if(in_iso_rate==0. || in_iso_rate>DBL_MAX)
    return in_iso_rate;

  if(iso_period_type==IsoPeriodType::constant)
    return in_iso_rate;

  // use collision-based isotropization
  double detP;
  if(iso_period_type==IsoPeriodType::det)
  {
    detP = det3x3symm(P11, P12, P13, P22, P23, P33);
  }
  else
  {
    assert_eq(iso_period_type,IsoPeriodType::trace);
    // use the determinant of the isotropized pressure
    detP = p*p*p;
  }
  // change to set an error flag?
  assert_gt(detP,0.);
  //if(detP<=0.) return  0.;
  const double rho5 = (rho*rho)*(rho*rho)*rho;
  const double pcl_mass_inv_3 = pcl_mass_inv*pcl_mass_inv*pcl_mass_inv;
  return sqrt(rho5/detP)*pcl_mass_inv_3*in_iso_rate;
}

// older version that works with periods instead of rates
//
inline double get_iso_period(
  IsoPeriodType::Enum iso_period_type, double in_iso_period, double pcl_mass,
  double rho, double p,
  double P11, double P12, double P13, double P22, double P23, double P33)
{
  if(!in_iso_period) return 0.;

  if(iso_period_type==IsoPeriodType::constant)
    return in_iso_period;

  // use collision-based isotropization
  double detP;
  if(iso_period_type==IsoPeriodType::det)
  {
    detP = det3x3symm(P11, P12, P13, P22, P23, P33);
  }
  else
  {
    assert_eq(iso_period_type,IsoPeriodType::trace);
    // use the determinant of the isotropized pressure
    detP = p*p*p;
  }
  // change to set an error flag?
  assert_gt(detP,0.);
  //if(detP<=0.) return  0.;
  const double rho5 = (rho*rho)*(rho*rho)*rho;
  const double pcl_mass3 = pcl_mass*pcl_mass*pcl_mass;
  return sqrt(detP/rho5)*pcl_mass3*in_iso_period;
}

inline void isotropizePressure(
    const double isotropization_rate,
    const double p_i,
    const double P11_i,
    const double P12_i,
    const double P13_i,
    const double P22_i,
    const double P23_i,
    const double P33_i,
    double& dP11_i,
    double& dP12_i,
    double& dP13_i,
    double& dP22_i,
    double& dP23_i,
    double& dP33_i)
{
    // const double p_i = (P11_i + P22_i + P33_i)*(1./3.0);
    dP11_i -= (P11_i-p_i)*isotropization_rate;
    dP12_i -= (P12_i    )*isotropization_rate;
    dP13_i -= (P13_i    )*isotropization_rate;
    dP22_i -= (P22_i-p_i)*isotropization_rate;
    dP23_i -= (P23_i    )*isotropization_rate;
    dP33_i -= (P33_i-p_i)*isotropization_rate;
}

inline void laterallyIsotropizePressure(
    const double isotropization_rate,
    const double B1,
    const double B2,
    const double B3,
    const double B,
    const double P11_i,
    const double P12_i,
    const double P13_i,
    const double P22_i,
    const double P23_i,
    const double P33_i,
    double& dP11_i,
    double& dP12_i,
    double& dP13_i,
    double& dP22_i,
    double& dP23_i,
    double& dP33_i)
{
    assert_gt(B,0.);
    // normalize the magnetic field
    // Maybe we should only partially relax
    // when the magnetic field is strong
    const double b1 = B1/B;
    const double b2 = B2/B;
    const double b3 = B2/B;
    const double bb11 = b1*b1;
    const double bb22 = b2*b2;
    const double bb33 = b3*b3;
    const double bb12 = b1*b2;
    const double bb13 = b1*b3;
    const double bb23 = b2*b3;

    // get pressure components
    //
    // parallel pressure component
    // p_parallel = b.P.b;
    const double p_parallel
      = bb11*P11_i
      + bb22*P22_i
      + bb33*P33_i
      + 2.*bb12*P12_i
      + 2.*bb13*P13_i
      + 2.*bb23*P23_i;
    // perpendicular pressure component
    // p_perp = P:(I_minus_bb)/2;
    const double p_perp = (
         (1.-bb11)*P11_i
       + (1.-bb22)*P22_i
       + (1.-bb33)*P33_i
       - 2.*P12_i*bb12
       - 2.*P13_i*bb13
       - 2.*P23_i*bb23);
    const double p_diff = p_parallel - p_perp;
    // P = p_parallel*bb + p_perp*(I-bb)
    //   = bb*(p_parallel-p_perp) + I*p_perp;
    const double lat_P11 = bb11*(p_diff) + p_perp ;
    const double lat_P22 = bb22*(p_diff) + p_perp ;
    const double lat_P33 = bb33*(p_diff) + p_perp ;
    const double lat_P12 = bb12*(p_diff) ;
    const double lat_P13 = bb13*(p_diff) ;
    const double lat_P23 = bb23*(p_diff) ;
    dP11_i -= (P11_i-lat_P11)*isotropization_rate;
    dP12_i -= (P12_i-lat_P12)*isotropization_rate;
    dP13_i -= (P13_i-lat_P13)*isotropization_rate;
    dP22_i -= (P22_i-lat_P22)*isotropization_rate;
    dP23_i -= (P23_i-lat_P23)*isotropization_rate;
    dP33_i -= (P33_i-lat_P33)*isotropization_rate;
}

// SourceSolver::relax_gas is a more up-to-date version of this. -eaj
//
// relax species s toward an isotropic pressure tensor
//
// note: q_new and qvals may actually point to the same data
inline bool isotropize_spc(
    dTensor2& q_new,
    const dTensor2& qvals,
    int numpts,
    double d_time,
    IsoPeriodType::Enum iso_period_type,
    double in_iso_period,
    double pcl_mass,
    int _rho_s,
    int _M1_s,
    int _M2_s,
    int _M3_s,
    int _N11_s,
    int _N12_s,
    int _N13_s,
    int _N22_s,
    int _N23_s,
    int _N33_s,
    int _B1_s,
    int _B2_s,
    int _B3_s)
{
    // determine whether/how to isotropize
    //
    const bool isotropize = (in_iso_period >= 0);
    if(!isotropize)
        return false;

    const bool physically_isotropize =
       (iso_period_type==IsoPeriodType::det
      ||iso_period_type==IsoPeriodType::trace);

    // instantaneously relax toward isotropy
    //
    // Loop over each point
    for (int i=1; i<=numpts; i++)
    {
        const double& rho = qvals.get(i, _rho_s);
        // change this: should not make a numerical assert
        // statement; instead should set an error flag which
        // is later checked and results in an exit with
        // proper behavior (e.g. dumping current state)
        // (probably easiest would be to throw and catch an error)
        assert_gt(rho,0.);
        const double& M1  = qvals.get(i, _M1_s);
        const double& M2  = qvals.get(i, _M2_s);
        const double& M3  = qvals.get(i, _M3_s);
        const double  u1  = M1/rho;
        const double  u2  = M2/rho;
        const double  u3  = M3/rho;
        const double& N11 = qvals.get(i, _N11_s);
        const double& N22 = qvals.get(i, _N22_s);
        const double& N33 = qvals.get(i, _N33_s);
        //
        const double D11 = M1*u1;
        const double D12 = M1*u2;
        const double D13 = M1*u3;
        const double D22 = M2*u2;
        const double D23 = M2*u3;
        const double D33 = M3*u3;
        //
        double P11 = N11-D11;
        double P22 = N22-D22;
        double P33 = N33-D33;
        double p = (P11 + P22 + P33)*(1./3.);
        //
        const double& N12 = qvals.get(i, _N12_s);
        const double& N13 = qvals.get(i, _N13_s);
        const double& N23 = qvals.get(i, _N23_s);
        //
        double P12 = N12-D12;
        double P13 = N13-D13;
        double P23 = N23-D23;

        // relax to isotropic pressure
        {
            double iso_period = Relax::get_iso_period(
              iso_period_type, in_iso_period, pcl_mass,
              rho, p, P11, P12, P13, P22, P23, P33);
            double portion_to_isotropize;
            if(iso_period > 0.)
              portion_to_isotropize = 1.-exp(-d_time/iso_period);
            else
              portion_to_isotropize = 1.;
            isotropizePressure(portion_to_isotropize, p,
                 P11, P12, P13, P22, P23, P33,
                 P11, P12, P13, P22, P23, P33);
            const bool laterally_isotropize = false;
            if(in_iso_period && laterally_isotropize)
            {
                const double& B1 = qvals.get(i,_B1_s);
                const double& B2 = qvals.get(i,_B2_s);
                const double& B3 = _B3_s>0 ? qvals.get(i,_B3_s) : 0.;
                const double B = sqrt(B1*B1 + B2*B2 + B3*B3);
                const double gyrofreq = B/pcl_mass;
                // const double base_lat_iso_timescale
                //   = plasmaParams.get_base_lat_iso_timescale();
                const double base_lat_iso_timescale = 1.;
                // it is not clear to me what lat_iso_rate should be
                const double lat_iso_rate = gyrofreq*
                  (base_lat_iso_timescale/iso_period);
                const double portion_to_gyrotropize = 1.-exp(-d_time*lat_iso_rate);
                
                if(portion_to_gyrotropize > 0.)
                {
                    Relax::laterallyIsotropizePressure(
                         portion_to_gyrotropize,
                         B1,B2,B3,B,
                         P11, P12, P13, P22, P23, P33,
                         P11, P12, P13, P22, P23, P33);
                }
            }
        }
        //
        q_new.set(i, _N11_s, D11 + P11);
        q_new.set(i, _N12_s, D12 + P12);
        q_new.set(i, _N13_s, D13 + P13);
        q_new.set(i, _N22_s, D22 + P22);
        q_new.set(i, _N23_s, D23 + P23);
        q_new.set(i, _N33_s, D33 + P33);
    }
    return true;
}

}
#endif
