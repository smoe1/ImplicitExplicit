#include "Components.h"
#include "tensors.h"
// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
void ProjectRightEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
    void ProjectRightEig_FiveMoment(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightEig_FiveMoment(ixy, 0, Q_ave, W, Q);

    void SymPair_ProjectRightEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    SymPair_ProjectRightEig_Maxwell(ixy, 5, Q_ave, W, Q);

    if(Q.getsize(2)<_y_i) return;
    void ProjectRightConvectedScalar(int idx,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightConvectedScalar(_x_i, Q_ave, W, Q);
    ProjectRightConvectedScalar(_y_i, Q_ave, W, Q);
}
