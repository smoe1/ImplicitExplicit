#include "SLState.h"
///////////////////////////////////////////////////////////////////////////////
// Function definitions for DogSolveSL.cpp
void CopyQ(const dTensorBC4&,dTensorBC4&);
double Max(const dTensorBC2& A);
void ConSoln(const dTensorBC4& aux, const dTensorBC4& q, double t);
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q);
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, void* data);
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q);
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q, void* data);
double GetSLCFL(double dt, const dTensorBC4& aux, const dTensorBC3& smax);
void StepAdvec(const double& dt,
        dTensorBC4& qold, dTensorBC4& qnew, 
        dTensorBC4& aux, const dTensorBC2& u1,
        const dTensorBC2& u2, int direction, SL_state& sl_state);
void StepAdvecHybrid(const double& dt,
        dTensorBC4& qold, dTensorBC4& qnew, 
        dTensorBC4& aux, const dTensorBC2& u1,
        const dTensorBC2& u2, int direction, SL_state& sl_state);
void SetAdvecSpeed(const dTensor2& phi, const dTensorBC4& q,
		   const dTensorBC4& aux, dTensorBC3& smax, const int &vel_dir,
		   dTensorBC2& u1, dTensorBC2& u2, SL_state& sl_state);
void SetAdvecSpeed_MOD(const dTensor2& phi1d, 
		       const dTensorBC4& q,
		       const dTensorBC3* Efield,
		       const dTensorBC3* v1d,
		       dTensorBC3& smax, 
		       const int &vel_dir,
		       dTensorBC2& u1, 
		       dTensorBC2& u2);
void ConstructL(
        const edge_data& EdgeData,
        dTensorBC4& aux, // SetBndValues modifies ghost cells
        dTensorBC4& q,   // SetBndValues modifies ghost cells
        dTensorBC4& Lstar, dTensorBC3& smax);
void ResInt(double dt, dTensorBC4* L[], dTensorBC5& ILout);
void TimeStepSDC(int morder, double t, double dt, 
        dTensor1& dtvec, dTensor1& tvec);
void ApplyScalarLimiter(const dTensorBC4& aux, dTensorBC4& q);
void ProjectRightEig(int,const dTensor1&,const dTensor1&,
        const dTensor2&,dTensor2&);
void ProjectLeftEig(int,const dTensor1&,const dTensor1&,
        const dTensor2&,dTensor2&);
void AddViscosity(const double& dt, const int& direction, 
    dTensorBC4& aux, dTensorBC4& q, dTensorBC4& Lstar);

// source term additions ... 
void SourceTermFunc_extra(const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& source, void *data);
void L2Project_extra(int istart, int iend, int jstart, int jend,
        const dTensorBC4& qin, const dTensorBC4& auxin, dTensorBC4& Fout,
        void (*SourceTermFunc)(const dTensor2& xpts,
            const dTensor2& qvals, const dTensor2& auxvals,
            dTensor2& source, void *data),
            void *data);
///////////////////////////////////////////////////////////////////////////////
