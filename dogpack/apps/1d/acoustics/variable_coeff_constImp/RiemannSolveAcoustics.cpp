#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"


//
double RiemannSolve(const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr,
        dTensor1& Fl,
        dTensor1& Fr,
        void (*FluxFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
        void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&))
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();

     double delta[2];
        delta[0] = Qr.get(1) - Ql.get(1);
        delta[1] = Qr.get(2) - Ql.get(2);


    //        # impedances:
    double   zr = Auxr.get(1)*Auxr.get(2);
    double   zl = Auxl.get(1)*Auxl.get(2);

    double    a1 = (-delta[0] + zr*delta[1]) / (zl + zr);
    double    a2 =  (delta[0] + zl*delta[1]) / (zl + zr);
    
    //        # Compute the waves.


    double Qm[2];
    Qm[0]=Ql.get(1)-a1*zl;
    Qm[1]=Ql.get(2)+a1;

    Fl.set(1, Auxl.get(1)*Auxl.get(1)*Auxl.get(2)*Qm[1]);
    Fl.set(2, Qm[0]/Auxl.get(2));

    Fr.set(1, Auxr.get(1)*Auxr.get(1)*Auxr.get(2)*Qm[1]);
    Fr.set(2, Qm[0]/Auxr.get(2));

    // Calculate minimum and maximum speeds
    double s1,s2;
    SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);
    

    double smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;
}
