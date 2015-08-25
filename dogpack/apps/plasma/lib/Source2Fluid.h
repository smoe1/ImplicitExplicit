#include "tensors.h"
#include "PlasmaParams.h"

namespace SourceSolver{

  // need to test this.
  void solve_electro_momentum_sympair(
    dTensor2& prim,
    const double d_time,
    int _rho_i,
    int _B1, int _B2, int _E3);

  void solve_electro_momentum(
    dTensor2& prim,
    const double d_time,
    int _rho_i, int _rho_e,
    int _B1, int _B2, int _B3,
    int _E1, int _E2, int _E3);

  void rotate_pressure_tensor(
    dTensor2& prim,
    double d_time,
    double q_over_m,
    int moffset,
    const int _B1_s,
    const int _B2_s,
    const int _B3_s);

  bool relax_gas(
    dTensor2& prim,
    double d_time,
    GasModelType::Enum gasModelType,
    SpeciesType::Enum speciesType,
    int moffset,
    const int _B1_s,
    const int _B2_s,
    const int _B3_s);

  // deprecated; calls relax_gas
  bool isotropize_pressure(
    dTensor2& prim,
    double d_time,
    SpeciesType::Enum speciesType,
    int moffset,
    const int _B1_s,
    const int _B2_s,
    const int _B3_s);
}
