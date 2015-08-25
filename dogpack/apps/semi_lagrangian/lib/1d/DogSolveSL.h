#ifndef __DOG_SOLVESL_H_
#define __DOG_SOLVESL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "tensors.h"
#include "DogParams.h"
using namespace std;

void DogSolveSL(const dTensor2& node, const dTensor1& prim_vol, dTensorBC3& aux, 
    dTensorBC3& qold, dTensorBC3& qnew, dTensorBC1& smax, double tstart, 
    double tend, int nv, const int method[], double dtv[], 
    const double cflv[], string outputdir);

void CopyQ(const dTensorBC3&,dTensorBC3&);
void BeforeStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
void BeforeStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&,void*);
void AfterStep(double,const dTensor2&,dTensorBC3&,dTensorBC3&);
double GetCFL(double,double,const dTensor1&,const int[],
        const dTensorBC3&,const dTensorBC1&);
void ConSoln(const int[],const dTensor2&,const dTensorBC3&,
        const dTensorBC3&,double,string);
// ------------------------------------------------------------

void StepAdvec(const double& dt, const dTensor2& node,
        dTensorBC3& auxvals, dTensorBC1& smax,
        dTensorBC3& qold,  // limiter will limit qold
        dTensorBC3& qnew);
void StepAdvecNonCons(const double& dt, const dTensor2& node,
        dTensorBC3& auxvals, dTensorBC1& smax,
        const dTensorBC3& qvals, dTensorBC3& qnew);
void StepAdvecFluxForm(const double& dt, const dTensor2& node,
        dTensorBC3& auxvals, dTensorBC1& smax,
        const dTensorBC3& qvals, dTensorBC3& qnew);

#endif
