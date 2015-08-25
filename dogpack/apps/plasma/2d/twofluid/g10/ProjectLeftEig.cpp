#include "Components.h"

// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
class dTensor1;
class dTensor2;
#include "tensors.h"
#include "assert.h"
void ProjectLeftEig(int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
  int n_offset_ions = 0;
  int n_offset_elcs = 10;
  int n_offset_maxwell = 20;

  // Choice of eigenstructure needs to match with choice in ProjectRightEig

  //void ProjectLeftEig_gas10(int ixy, int moffset,
  //  const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  //ProjectLeftEig_gas10(ixy, n_offset_ions, Q_ave, Q, W);
  //ProjectLeftEig_gas10(ixy, n_offset_elcs, Q_ave, Q, W);

  void ProjectLeftEig_10moment(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& WJ);
  ProjectLeftEig_10moment(ixy, n_offset_ions, Q_ave, Q, W);
  ProjectLeftEig_10moment(ixy, n_offset_elcs, Q_ave, Q, W);

  // check my eigenvectors against James
  #if 0
  {
    // Evidently 5 and 6 to not match with James's 6 and 5.
    int A2J[11] = {0, 1, 2, 3, 4, 6, 5, 7, 9, 8, 10};
    dTensor2 WJ(W.getsize(1),W.getsize(2));
    void ProjectLeftEig_10moment(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& WJ);
    ProjectLeftEig_10moment(ixy, n_offset_ions, Q_ave, Q, WJ);
    ProjectLeftEig_10moment(ixy, n_offset_elcs, Q_ave, Q, WJ);
    for(int k=1;k<=W.getsize(2);k++)
    for(int m=1;m<=10;m++)
    {
      if(fcmp(WJ.get(A2J[m]+n_offset_ions,k),W.get(m+n_offset_ions,k),1e-13))
      {
        printf("WJ.get(%d+n_offset_ions,%d)=%24.16e\n",m,k,WJ.get(A2J[m]+n_offset_ions,k));
        printf(" W.get(%d+n_offset_ions,%d)=%24.16e\n",m,k, W.get(m+n_offset_ions,k));
      }
      //assert_almost_eq(WJ.get(m+n_offset_ions,k),W.get(m+n_offset_ions,k));
      //assert_almost_eq(WJ.get(m+n_offset_elcs,k),W.get(m+n_offset_elcs,k));
    }
  }
  #endif

  void ProjectLeftEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftEig_Maxwell(ixy, n_offset_maxwell, Q_ave, Q, W);

  void ProjectLeftConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
  ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
  ProjectLeftConvectedScalar(_entropy_e, Q_ave, Q, W);
}
