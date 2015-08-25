#include "tensors.h"
#include "Params.h"

// instantaneously relax
//
void Relax(const dTensor1& xpts, 
	   const dTensor2& qvals, 
	   const dTensor2& auxvals, 
	   dTensor2& q_new)
{
    int i,k;
    int numpts=xpts.getsize();

    // Loop over each point
    for (i=1; i<=numpts; i++)
    {
        // relax ions
        //
        // really ought to test whether
        // the relaxation is too stiff and if so
        // relax instantaneously rather than
        // using the source term.
        //
        if(appParams.ion_rlx_period==0.0)
        {
            const double& rho_i = qvals.get(i,1);
            const double& M1_i  = qvals.get(i,2);
            const double& M2_i  = qvals.get(i,3);
            const double& M3_i  = qvals.get(i,4);
            const double  u1_i  = M1_i/rho_i;
            const double  u2_i  = M2_i/rho_i;
            const double  u3_i  = M3_i/rho_i;
            const double& N11_i = qvals.get(i,5);
            const double& N12_i = qvals.get(i,6);
            const double& N13_i = qvals.get(i,7);
            const double& N22_i = qvals.get(i,8);
            const double& N23_i = qvals.get(i,9);
            const double& N33_i = qvals.get(i,10);
            //
            const double D11_i = M1_i*u1_i;
            const double D12_i = M1_i*u2_i;
            const double D13_i = M1_i*u3_i;
            const double D22_i = M2_i*u2_i;
            const double D23_i = M2_i*u3_i;
            const double D33_i = M3_i*u3_i;
            //
            double P11_i = N11_i-D11_i;
            double P12_i = N12_i-D12_i;
            double P13_i = N13_i-D13_i;
            double P22_i = N22_i-D22_i;
            double P23_i = N23_i-D23_i;
            double P33_i = N33_i-D33_i;
            double p_i   = (P11_i + P22_i + P33_i)/3.0;
            //
            // Instantaneously relax to isotropic pressure
            //
            P11_i = p_i;
            P12_i = 0.0;
            P13_i = 0.0;
            P22_i = p_i;
            P23_i = 0.0;
            P33_i = p_i;
            //
            // this copying actually does nothing,
            // because qvals and q_new reference the same memory.
            for(k=1;k<5;k++) q_new.set(i,k, qvals.get(i,k));
            //
            q_new.set(i,5,  P11_i + D11_i); // N11_i 
            q_new.set(i,6,  P12_i + D12_i); // N12_i 
            q_new.set(i,7,  P13_i + D13_i); // N13_i 
            q_new.set(i,8,  P22_i + D22_i); // N22_i 
            q_new.set(i,9,  P23_i + D23_i); // N23_i 
            q_new.set(i,10, P33_i + D33_i); // N33_i 
        }
        else
        {
            // this copying actually does nothing
            for(k=1;k<11;k++) q_new.set(i,k, qvals.get(i,k));
        }
        
        if(appParams.elc_rlx_period==0.0)
        {
            const double& rho_e = qvals.get(i,11);
            const double& M1_e  = qvals.get(i,12);
            const double& M2_e  = qvals.get(i,13);
            const double& M3_e  = qvals.get(i,14);
            const double  u1_e  = M1_e/rho_e;
            const double  u2_e  = M2_e/rho_e;
            const double  u3_e  = M3_e/rho_e;
            const double& N11_e = qvals.get(i,15);
            const double& N12_e = qvals.get(i,16);
            const double& N13_e = qvals.get(i,17);
            const double& N22_e = qvals.get(i,18);
            const double& N23_e = qvals.get(i,19);
            const double& N33_e = qvals.get(i,20);
            //
            const double D11_e = M1_e*u1_e;
            const double D12_e = M1_e*u2_e;
            const double D13_e = M1_e*u3_e;
            const double D22_e = M2_e*u2_e;
            const double D23_e = M2_e*u3_e;
            const double D33_e = M3_e*u3_e;
            //
            double P11_e = N11_e-D11_e;
            double P12_e; // = N12_e-D12_e;
            double P13_e; // = N13_e-D13_e;
            double P22_e = N22_e-D22_e;
            double P23_e; // = N23_e-D23_e;
            double P33_e = N33_e-D33_e;
            double p_e   = (P11_e + P22_e + P33_e)/3.0;
            //
            // Instantaneously relax to isotropic pressure
            //
            P11_e = p_e;
            P12_e = 0.0;
            P13_e = 0.0;
            P22_e = p_e;
            P23_e = 0.0;
            P33_e = p_e;
            //
            // this copying actually does nothing
            for(k=11;k<15;k++) q_new.set(i,k, qvals.get(i,k));
            //
            q_new.set(i,15,  P11_e + D11_e); // N11_e 
            q_new.set(i,16,  P12_e + D12_e); // N12_e 
            q_new.set(i,17,  P13_e + D13_e); // N13_e 
            q_new.set(i,18,  P22_e + D22_e); // N22_e 
            q_new.set(i,19,  P23_e + D23_e); // N23_e 
            q_new.set(i,20,  P33_e + D33_e); // N33_e 
        }
        else // do not relax electrons.
        {
            // this copying actually does nothing
            for(k=11;k<27;k++) q_new.set(i,k, qvals.get(i,k));
        }
        
        // this copying actually does nothing
        for(k=21;k<27;k++) q_new.set(i,k, qvals.get(i,k));
    }
}

void AfterStep(double dt, 
	       const dTensor2& node, 
	       dTensorBC3& aux, 
	       dTensorBC3& q)
{
    int i,j;
    int melems = q.getsize(1);
    int   meqn = q.getsize(2);
    int   kmax = q.getsize(3);
    int   maux = aux.getsize(2);
    int    mbc = q.getmbc();
    void L2Project(int mopt, int istart, int iend,
		   const dTensor2& node,
		   const dTensorBC3& qin, 
		   const dTensorBC3& auxin,  
		   dTensorBC3& Fout,
		   void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
   
    L2Project(0,1-mbc,melems+mbc,node,q,aux,q,&Relax);
}
