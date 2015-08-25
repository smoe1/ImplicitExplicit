#include <cmath>
#include "constants.h"
#include "tensors.h"
#include "Params.h"

// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& W,
		     dTensor2& Q)
{    
    int m,k;
    int meqn = Q.getsize(1);
    int kmax = Q.getsize(2)+1;
    double rho,u1,u2,u3,P11,P12,P13,P22,P23,P33,c;
    dTensor2 Rm(10,10);

    double& cs_light = appParams.cs_light;
  
    // -----------------------------
    // For the IONS
    // -----------------------------

    // Average states
    rho = Q_ave.get(1);
    u1  = Q_ave.get(2)/rho;
    u2  = Q_ave.get(3)/rho;
    u3  = Q_ave.get(4)/rho;
    P11 = Q_ave.get(5)  - rho*u1*u1;
    P12 = Q_ave.get(6)  - rho*u1*u2;
    P13 = Q_ave.get(7)  - rho*u1*u3;
    P22 = Q_ave.get(8)  - rho*u2*u2;
    P23 = Q_ave.get(9)  - rho*u2*u3;
    P33 = Q_ave.get(10) - rho*u3*u3;
    c   = sqrt(P11/rho);

    // Set non-zero elements of matrix of right eigenvectors
    Rm.set(1,1,        P11*rho );
    Rm.set(2,1,     u1*P11*rho - rho*sq3*c*P11 );
    Rm.set(3,1,     u2*P11*rho - rho*sq3*c*P12 );
    Rm.set(4,1,     u3*P11*rho - rho*sq3*c*P13 );
    Rm.set(5,1,  u1*u1*P11*rho - 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
    Rm.set(6,1,  u1*u2*P11*rho - rho*u2*sq3*c*P11 - rho*u1*sq3*c*P12 + 3.0*P11*P12 );
    Rm.set(7,1,  u1*u3*P11*rho - rho*u3*sq3*c*P11 - rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
    Rm.set(8,1,  u2*u2*P11*rho - 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
    Rm.set(9,1,  u2*u3*P11*rho - rho*u3*sq3*c*P12 - rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
    Rm.set(10,1, u3*u3*P11*rho - 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
    
    Rm.set(3,2,  -rho*c );
    Rm.set(6,2,  -rho*u1*c + P11 );
    Rm.set(8,2,  -2.0*rho*u2*c + 2.0*P12 );
    Rm.set(9,2,  -rho*u3*c + P13 );

    Rm.set(4,3,  -rho*c ); 
    Rm.set(7,3,  -rho*u1*c + P11 );
    Rm.set(9,3,  -rho*u2*c + P12 );
    Rm.set(10,3, -2.0*rho*u3*c + 2.0*P13 );
    
    Rm.set(1,4,  1.0 );
    Rm.set(2,4,  u1 );
    Rm.set(3,4,  u2 );
    Rm.set(4,4,  u3 );
    Rm.set(5,4,  u1*u1 );
    Rm.set(6,4,  u1*u2 ); 
    Rm.set(7,4,  u1*u3 ); 
    Rm.set(8,4,  u2*u2 );
    Rm.set(9,4,  u2*u3 ); 
    Rm.set(10,4, u3*u3 );
    
    Rm.set(8,5,  1.0 );
    
    Rm.set(9,6,  1.0 );     

    Rm.set(10,7, 1.0 );
    
    Rm.set(3,8,  rho*c );
    Rm.set(6,8,  rho*u1*c + P11 );
    Rm.set(8,8,  2.0*rho*u2*c + 2.0*P12 );
    Rm.set(9,8,  rho*u3*c + P13 );
    
    Rm.set(4,9,  rho*c ); 
    Rm.set(7,9,  rho*u1*c + P11 );
    Rm.set(9,9,  rho*u2*c + P12 );
    Rm.set(10,9, 2.0*rho*u3*c + 2.0*P13 );

    Rm.set(1,10,        P11*rho );
    Rm.set(2,10,     u1*P11*rho + rho*sq3*c*P11 );
    Rm.set(3,10,     u2*P11*rho + rho*sq3*c*P12 );
    Rm.set(4,10,     u3*P11*rho + rho*sq3*c*P13 );
    Rm.set(5,10,  u1*u1*P11*rho + 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
    Rm.set(6,10,  u1*u2*P11*rho + rho*u2*sq3*c*P11 + rho*u1*sq3*c*P12 + 3.0*P11*P12 );
    Rm.set(7,10,  u1*u3*P11*rho + rho*u3*sq3*c*P11 + rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
    Rm.set(8,10,  u2*u2*P11*rho + 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
    Rm.set(9,10,  u2*u3*P11*rho + rho*u3*sq3*c*P12 + rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
    Rm.set(10,10, u3*u3*P11*rho + 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Q.set(1,k, Rm.get(1,1)*W.get(1,k) + W.get(4,k) + Rm.get(1,10)*W.get(10,k) );

        Q.set(2,k, Rm.get(2,1)*W.get(1,k) + Rm.get(2,4)*W.get(4,k) + Rm.get(2,10)*W.get(10,k) );

	Q.set(3,k, Rm.get(3,1)*W.get(1,k) + Rm.get(3,2)*W.get(2,k) + Rm.get(3,4)*W.get(4,k) 
	      + Rm.get(3,8)*W.get(8,k) + Rm.get(3,10)*W.get(10,k) );
	
	Q.set(4,k, Rm.get(4,1)*W.get(1,k) + Rm.get(4,3)*W.get(3,k) + Rm.get(4,4)*W.get(4,k) 
	      + Rm.get(4,9)*W.get(9,k) + Rm.get(4,10)*W.get(10,k) );

	Q.set(5,k, Rm.get(5,1)*W.get(1,k) + Rm.get(5,4)*W.get(4,k) + Rm.get(5,10)*W.get(10,k) );

	Q.set(6,k, Rm.get(6,1)*W.get(1,k) + Rm.get(6,2)*W.get(2,k) + Rm.get(6,4)*W.get(4,k) 
	      + Rm.get(6,8)*W.get(8,k) + Rm.get(6,10)*W.get(10,k) );

	Q.set(7,k, Rm.get(7,1)*W.get(1,k) + Rm.get(7,3)*W.get(3,k) + Rm.get(7,4)*W.get(4,k) 
	      + Rm.get(7,9)*W.get(9,k) + Rm.get(7,10)*W.get(10,k) );

	Q.set(8,k, Rm.get(8,1)*W.get(1,k) + Rm.get(8,2)*W.get(2,k) + Rm.get(8,4)*W.get(4,k) 
	      + W.get(5,k) + Rm.get(8,8)*W.get(8,k)  + Rm.get(8,10)*W.get(10,k) );

	Q.set(9,k, Rm.get(9,1)*W.get(1,k) + Rm.get(9,2)*W.get(2,k) + Rm.get(9,3)*W.get(3,k) 
	      + Rm.get(9,4)*W.get(4,k) + W.get(6,k) + Rm.get(9,8)*W.get(8,k) 
	      + Rm.get(9,9)*W.get(9,k) + Rm.get(9,10)*W.get(10,k) );

	Q.set(10,k, Rm.get(10,1)*W.get(1,k) + Rm.get(10,3)*W.get(3,k) + Rm.get(10,4)*W.get(4,k) 
	      + W.get(7,k) + Rm.get(10,9)*W.get(9,k) + Rm.get(10,10)*W.get(10,k) );
    }


    // -----------------------------
    // For the ELECTRONS
    // -----------------------------

    // Average states
    rho = Q_ave.get(11);
    u1  = Q_ave.get(12)/rho;
    u2  = Q_ave.get(13)/rho;
    u3  = Q_ave.get(14)/rho;
    P11 = Q_ave.get(15) - rho*u1*u1;
    P12 = Q_ave.get(16) - rho*u1*u2;
    P13 = Q_ave.get(17) - rho*u1*u3;
    P22 = Q_ave.get(18) - rho*u2*u2;
    P23 = Q_ave.get(19) - rho*u2*u3;
    P33 = Q_ave.get(20) - rho*u3*u3;
    c   = sqrt(P11/rho);

    // Set non-zero elements of matrix of right eigenvectors
    Rm.set(1,1,        P11*rho );
    Rm.set(2,1,     u1*P11*rho - rho*sq3*c*P11 );
    Rm.set(3,1,     u2*P11*rho - rho*sq3*c*P12 );
    Rm.set(4,1,     u3*P11*rho - rho*sq3*c*P13 );
    Rm.set(5,1,  u1*u1*P11*rho - 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
    Rm.set(6,1,  u1*u2*P11*rho - rho*u2*sq3*c*P11 - rho*u1*sq3*c*P12 + 3.0*P11*P12 );
    Rm.set(7,1,  u1*u3*P11*rho - rho*u3*sq3*c*P11 - rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
    Rm.set(8,1,  u2*u2*P11*rho - 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
    Rm.set(9,1,  u2*u3*P11*rho - rho*u3*sq3*c*P12 - rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
    Rm.set(10,1, u3*u3*P11*rho - 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
    
    Rm.set(3,2,  -rho*c );
    Rm.set(6,2,  -rho*u1*c + P11 );
    Rm.set(8,2,  -2.0*rho*u2*c + 2.0*P12 );
    Rm.set(9,2,  -rho*u3*c + P13 );

    Rm.set(4,3,  -rho*c ); 
    Rm.set(7,3,  -rho*u1*c + P11 );
    Rm.set(9,3,  -rho*u2*c + P12 );
    Rm.set(10,3, -2.0*rho*u3*c + 2.0*P13 );
    
    Rm.set(1,4,  1.0 );
    Rm.set(2,4,  u1 );
    Rm.set(3,4,  u2 );
    Rm.set(4,4,  u3 );
    Rm.set(5,4,  u1*u1 );
    Rm.set(6,4,  u1*u2 ); 
    Rm.set(7,4,  u1*u3 ); 
    Rm.set(8,4,  u2*u2 );
    Rm.set(9,4,  u2*u3 ); 
    Rm.set(10,4, u3*u3 );
    
    Rm.set(8,5,  1.0 );
    
    Rm.set(9,6,  1.0 );     

    Rm.set(10,7, 1.0 );
    
    Rm.set(3,8,  rho*c );
    Rm.set(6,8,  rho*u1*c + P11 );
    Rm.set(8,8,  2.0*rho*u2*c + 2.0*P12 );
    Rm.set(9,8,  rho*u3*c + P13 );
    
    Rm.set(4,9,  rho*c ); 
    Rm.set(7,9,  rho*u1*c + P11 );
    Rm.set(9,9,  rho*u2*c + P12 );
    Rm.set(10,9, 2.0*rho*u3*c + 2.0*P13 );

    Rm.set(1,10,        P11*rho );
    Rm.set(2,10,     u1*P11*rho + rho*sq3*c*P11 );
    Rm.set(3,10,     u2*P11*rho + rho*sq3*c*P12 );
    Rm.set(4,10,     u3*P11*rho + rho*sq3*c*P13 );
    Rm.set(5,10,  u1*u1*P11*rho + 2.0*rho*u1*sq3*c*P11 + 3.0*P11*P11 );
    Rm.set(6,10,  u1*u2*P11*rho + rho*u2*sq3*c*P11 + rho*u1*sq3*c*P12 + 3.0*P11*P12 );
    Rm.set(7,10,  u1*u3*P11*rho + rho*u3*sq3*c*P11 + rho*u1*sq3*c*P13 + 3.0*P11*P13 ); 
    Rm.set(8,10,  u2*u2*P11*rho + 2.0*rho*u2*sq3*c*P12 + P11*P22 + 2.0*P12*P12 );
    Rm.set(9,10,  u2*u3*P11*rho + rho*u3*sq3*c*P12 + rho*u2*sq3*c*P13 + P11*P23 + 2.0*P13*P12 );
    Rm.set(10,10, u3*u3*P11*rho + 2.0*rho*u3*sq3*c*P13 + P11*P33 + 2.0*P13*P13 );
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Q.set(11,k, Rm.get(1,1)*W.get(11,k) + W.get(14,k) + Rm.get(1,10)*W.get(20,k) );

        Q.set(12,k, Rm.get(2,1)*W.get(11,k) + Rm.get(2,4)*W.get(14,k) + Rm.get(2,10)*W.get(20,k) );

	Q.set(13,k, Rm.get(3,1)*W.get(11,k) + Rm.get(3,2)*W.get(12,k) + Rm.get(3,4)*W.get(14,k) 
	      + Rm.get(3,8)*W.get(18,k) + Rm.get(3,10)*W.get(20,k) );
	
	Q.set(14,k, Rm.get(4,1)*W.get(11,k) + Rm.get(4,3)*W.get(13,k) + Rm.get(4,4)*W.get(14,k) 
	      + Rm.get(4,9)*W.get(19,k) + Rm.get(4,10)*W.get(20,k) );

	Q.set(15,k, Rm.get(5,1)*W.get(11,k) + Rm.get(5,4)*W.get(14,k) + Rm.get(5,10)*W.get(20,k) );

	Q.set(16,k, Rm.get(6,1)*W.get(11,k) + Rm.get(6,2)*W.get(12,k) + Rm.get(6,4)*W.get(14,k) 
	      + Rm.get(6,8)*W.get(18,k) + Rm.get(6,10)*W.get(20,k) );

	Q.set(17,k, Rm.get(7,1)*W.get(11,k) + Rm.get(7,3)*W.get(13,k) + Rm.get(7,4)*W.get(14,k) 
	      + Rm.get(7,9)*W.get(19,k) + Rm.get(7,10)*W.get(20,k) );

	Q.set(18,k, Rm.get(8,1)*W.get(11,k) + Rm.get(8,2)*W.get(12,k) + Rm.get(8,4)*W.get(14,k) 
	      + W.get(15,k) + Rm.get(8,8)*W.get(18,k)  + Rm.get(8,10)*W.get(20,k) );

	Q.set(19,k, Rm.get(9,1)*W.get(11,k) + Rm.get(9,2)*W.get(12,k) + Rm.get(9,3)*W.get(13,k) 
	      + Rm.get(9,4)*W.get(14,k) + W.get(16,k) + Rm.get(9,8)*W.get(18,k) 
	      + Rm.get(9,9)*W.get(19,k) + Rm.get(9,10)*W.get(20,k) );

	Q.set(20,k, Rm.get(10,1)*W.get(11,k) + Rm.get(10,3)*W.get(13,k) + Rm.get(10,4)*W.get(14,k) 
	      + W.get(17,k) + Rm.get(10,9)*W.get(19,k) + Rm.get(10,10)*W.get(20,k) );
    }


    // -----------------------------
    // For the ELECTROMAGNETIC FIELD
    // -----------------------------
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Q.set(21,k, W.get(23,k) );

	Q.set(22,k, W.get(21,k) + W.get(25,k) );

	Q.set(23,k, W.get(22,k) + W.get(26,k) );
	
	Q.set(24,k, W.get(24,k) );

	Q.set(25,k, cs_light*( W.get(26,k) - W.get(22,k) ) );

	Q.set(26,k, cs_light*( W.get(21,k) - W.get(25,k) ) );
    }
}
