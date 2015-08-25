#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
  int n_offset_ions = 0;
  int n_offset_maxwell = 10;

  //void ProjectRightEig_10moment(int ixy, int n_offset,
  //    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  //ProjectRightEig_10moment(ixy, n_offset_ions, Q_ave, W, Q);

  // limiting gas dynamics in conserved variables
  for (int k=1; k<=Q.getsize(2); k++)
  for (int n=1; n<=10; n++)
  {
      Q.set(n,k, W.get(n,k));
  }

  void SymPair_ProjectRightEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  SymPair_ProjectRightEig_Maxwell(ixy, n_offset_maxwell, Q_ave, W, Q);

  void ProjectRightConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  if(Q.getsize(1)<_entropy_i) return;
  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
}

