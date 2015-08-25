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
  int n_offset_elcs = 10;
  int n_offset_maxwell = 20;

  // Choice of eigenstructure needs to match with choice in ProjectLeftEig

  // Alec's eigenstructure (gives magnetic islands)
  // 
  //void ProjectRightEig_gas10(int ixy, int n_offset,
  //    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  //ProjectRightEig_gas10(ixy, n_offset_ions, Q_ave, W, Q);
  //ProjectRightEig_gas10(ixy, n_offset_elcs, Q_ave, W, Q);

  // James's eigenstructure (slower)
  //
  void ProjectRightEig_10moment(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& QJ);
  ProjectRightEig_10moment(ixy, n_offset_ions, Q_ave, W, Q);
  ProjectRightEig_10moment(ixy, n_offset_elcs, Q_ave, W, Q);

  #if 0
  {
    // Evidently 5 and 6 to not match with James's 6 and 5.
    int A2J[11] = {0, 1, 2, 3, 4, 6, 5, 7, 9, 8, 10};
    dTensor2 QJ(Q.getsize(1),Q.getsize(2));

    void ProjectRightEig_10moment(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& QJ);
    ProjectRightEig_10moment(ixy, n_offset_ions, Q_ave, W, QJ);
    ProjectRightEig_10moment(ixy, n_offset_elcs, Q_ave, W, QJ);

    for(int k=1;k<=Q.getsize(2);k++)
    for(int m=1;m<=10;m++)
    {
      if(fcmp(QJ.get(m+n_offset_ions,k),Q.get(m+n_offset_ions,k),1e-13))
      {
        printf("QJ.get(%d+n_offset_ions,%d)=%24.16e\n",m,k,QJ.get(m+n_offset_ions,k));
        printf(" Q.get(%d+n_offset_ions,%d)=%24.16e\n",m,k, Q.get(m+n_offset_ions,k));
      }
      //assert_almost_eq(QJ.get(m+n_offset_ions,k),Q.get(m+n_offset_ions,k));
      //assert_almost_eq(QJ.get(m+n_offset_elcs,k),Q.get(m+n_offset_elcs,k));
    }
  }
  #endif

  void ProjectRightEig_Maxwell(int ixy, int n_offset,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightEig_Maxwell(ixy, n_offset_maxwell, Q_ave, W, Q);

  void ProjectRightConvectedScalar(int idx,
      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
  ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);
}

