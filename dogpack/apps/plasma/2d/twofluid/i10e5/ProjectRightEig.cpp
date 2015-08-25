#include "Components.h"
// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
class dTensor1;
class dTensor2;
void ProjectRightEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
  int n_offset_ions = 0;
  int n_offset_elcs = 10;
  int n_offset_maxwell = 15;

  void ProjectRightEig_10moment(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightEig_10moment(ixy, n_offset_ions, Q_ave, W, Q);
  void ProjectRightEig_FiveMoment(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightEig_FiveMoment(ixy, n_offset_elcs, Q_ave, W, Q);

  void ProjectRightEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightEig_Maxwell(ixy, n_offset_maxwell, Q_ave, W, Q);

  void ProjectRightConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
  ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);
}

