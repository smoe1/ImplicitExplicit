#include "Components.h"
// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
class dTensor1;
class dTensor2;
void ProjectLeftEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
  int n_offset_ions = 0;
  int n_offset_elcs = 10;
  int n_offset_maxwell = 15;

  void ProjectLeftEig_10moment(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftEig_10moment(ixy, n_offset_ions, Q_ave, Q, W);
  void ProjectLeftEig_FiveMoment( int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftEig_FiveMoment(ixy, n_offset_elcs, Q_ave, Q, W);

  void ProjectLeftEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftEig_Maxwell(ixy, n_offset_maxwell, Q_ave, Q, W);

  void ProjectLeftConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
  ProjectLeftConvectedScalar(_entropy_e, Q_ave, Q, W);
}
