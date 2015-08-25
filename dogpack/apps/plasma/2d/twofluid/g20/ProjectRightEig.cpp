#include "Components.h"
// This is a user-supplied routine that projects
// W onto the right eigenvectors of the flux 
// Jacobian; the result is stored in Q
//
class dTensor1;
class dTensor2;
#include "tensors.h"
#include "assert.h"
void ProjectRightEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
  int n_offset_ions = 0;
  int n_offset_elcs = 20;
  int n_offset_maxwell = 40;

  //void ProjectRightEig_10moment(int ixy, int n_offset,
  //    const dTensor1& Q_ave, const dTensor2& W, dTensor2& QJ);
  //ProjectRightEig_10moment(ixy, n_offset_ions, Q_ave, W, Q);
  //ProjectRightEig_10moment(ixy, n_offset_elcs, Q_ave, W, Q);

  // limiting gas dynamics in conserved variables
  for (int k=1; k<=Q.getsize(2); k++)
  for (int n=1; n<=40; n++)
  {
      Q.set(n,k, W.get(n,k));
  }

  void ProjectRightEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightEig_Maxwell(ixy, n_offset_maxwell, Q_ave, W, Q);

  void ProjectRightConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
  ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);
}

