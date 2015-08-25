#ifndef _RKSOLVE1D_H_
#define _RKSOLVE1D_H_

#include <cmath>
#include <iostream>
#include <pthread.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "discontCell.h"
#include "SLState.h"
#include "tensors.h"
#include "dog_math.h"

//void CopyQ(const dTensorBC3& qin, dTensorBC3& qout );

double RKStep1D(
    double t,
    double dt, double xlow, double dx,
    const dTensor1& speeds,
    dTensorBC3& qold,  // Set Boundary conditions modifies qold ...
    dTensorBC3& qnew);

void ConstructL1D(
        double t, 
        double sgn_dt,
        const double xlow, const double dx,
        const dTensor1& speeds,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC3& Lstar );

double EulerStep( double t, double dt, 
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& qnew );

void SetPeriodicBndy1D( dTensorBC3 &qnew );

void L2ProjectAdvection(int mopt, int istart, int iend,
        double dx,
        const dTensorBC3& qin, 
        const dTensor1& speeds,  
        dTensorBC3& Fout
    );

void L2ProjectSource(double t, int istart, int iend,
        double dx, 
        const dTensor1& speeds,  
        dTensorBC3& Fout,
        void (*Func)(
            double t, 
            const dTensor1& speeds, 
            const dTensor1& xpts,
            dTensor2& fvals));

void HybridSourceTermFunc1D(double t, const dTensor1& speeds, const dTensor1&
    xpts, dTensor2& fvals);

#endif
