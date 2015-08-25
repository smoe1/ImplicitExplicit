#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
    void ProjectLeftEig_FiveMoment( int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftEig_FiveMoment(ixy, 0, Q_ave, Q, W);

    void SymPair_ProjectLeftEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    SymPair_ProjectLeftEig_Maxwell(ixy, 5, Q_ave, Q, W);

    if(Q.getsize(2)<_y_i) return;
    void ProjectLeftConvectedScalar(int idx,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftConvectedScalar(_x_i, Q_ave, Q, W);
    ProjectLeftConvectedScalar(_y_i, Q_ave, Q, W);
}
