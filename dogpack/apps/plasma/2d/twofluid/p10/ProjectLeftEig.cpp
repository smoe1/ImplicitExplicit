#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
void ProjectLeftEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
  int n_offset_ions = 0;
  int n_offset_maxwell = 10;

  //void ProjectLeftEig_10moment(int ixy, int n_offset,
  //    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  //ProjectLeftEig_10moment(ixy, n_offset_ions, Q_ave, Q, W);

  // limit gas dynamics in conserved variables
  for (int k=1; k<=W.getsize(2); k++)
  for (int n=1; n<=10; n++)
  {
      W.set(n,k, Q.get(n,k));
  }

  void SymPair_ProjectLeftEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  SymPair_ProjectLeftEig_Maxwell(ixy, n_offset_maxwell, Q_ave, Q, W);

  void ProjectLeftConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  if(Q.getsize(1)<_entropy_i) return;
  ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
}
