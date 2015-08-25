#ifndef __STEP_ADVECHYBRID_H_
#define __STEP_ADVECHYBRID_H_

#include <math.h>
#include <iostream>
#include <pthread.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "discontCell.h"
#include "SLState.h"

void StepAdvec1(int j_start, int j_end, const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew,  
        const dTensorBC2& u1, const dTensorBC2& u2, SL_state& sl_state);
void StepAdvec2(int i_start, int i_end, const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2 );

void StepAdvec1lines(int j_start, int j_end, const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2, SL_state& sl_state);

void StepAdvec2lines(int i_start, int i_end, const double& dt, 
        dTensorBC4& qold, dTensorBC4& qnew, 
        const dTensorBC2& u1, const dTensorBC2& u2);


// This routine is used for 1D advections 
void StepAdvec1D   ( double dt, double xlow, double dx, const dTensor1& speeds, const dTensorBC3& qold, dTensorBC3& qnew);
void RKSolve1D     ( double t, double dt, double xlow, double dx, const dTensor1& speeds, dTensorBC3& qold, dTensorBC3& qnew);
void CopyQ(const dTensorBC4& qin, dTensorBC4& qout );

#endif
