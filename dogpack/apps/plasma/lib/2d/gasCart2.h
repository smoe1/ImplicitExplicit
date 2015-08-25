#ifndef gasCart2_h
#define gasCart2_h

class dTensor2;
class dTensorBC4;
void UpdateEntropies(
    double dt,dTensorBC4& aux,
    dTensorBC4& q,
    int firstEntropyIdx,
    int numEntropyIndices,
    void (*reset_entropy)(
       const dTensor2& xpts,
       const dTensor2& qvals,
       const dTensor2& auxvals,
       dTensor2& q_new)
    );

#endif
