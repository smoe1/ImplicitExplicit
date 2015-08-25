#include <cmath>
#include "constants.h"
#include "tensors.h"
#include "Params.h"

// This is a user-supplied routine that projects
// Q onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in W
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
		    const dTensor1& Q_ave,
		    const dTensor2& Q,
		    dTensor2& W)
{    
    int m,k;
    int meqn = Q.getsize(1);
    int kmax = Q.getsize(2)+1;
    double rho,u1,u2,u3,P11,P12,P13,P22,P23,P33,c;
    dTensor2 Lm(10,10);

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

    // Set non-zero elements of matrix of left eigenvectors
    Lm.set(1,1,   u1*(sq3*P11 + rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,2,  -(sq3*P11 + 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,5,   1.0/(6.0*P11*P11) );
    
    Lm.set(2,1,  -(P11*u1*P12 - u2*P11*P11 + P12*u1*u1*c*rho - u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,2,   (P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,3,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(2,5,  -P12/(2.0*P11*P11) );
    Lm.set(2,6,   1.0/(2.0*P11) );

    Lm.set(3,1,  -(P11*P13*u1 - u3*P11*P11 + P13*u1*u1*c*rho - u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,2,   (P11*P13 + 2.0*u1*P13*c*rho - u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,4,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(3,5,  -P13/(2.0*P11*P11) );
    Lm.set(3,7,   1.0/(2.0*P11) );
    
    Lm.set(4,1,  -(rho*u1*u1 - 3.0*P11)/(3.0*P11) );
    Lm.set(4,2,   2.0*rho*u1/(3.0*P11) );
    Lm.set(4,5,  -rho/(3.0*P11) );

    Lm.set(5,1,  -(u1*u1*P11*P22 - 4.0*u1*u1*P12*P12 + 6.0*u1*P12*u2*P11 - 3.0*u2*u2*P11*P11)/(3.0*P11*P11) );
    Lm.set(5,2,   2.0*(P11*P22*u1 - 4.0*u1*P12*P12 + 3.0*u2*P11*P12)/(3.0*P11*P11) );
    Lm.set(5,3,   (2.0*u1*P12 - 2.0*u2*P11)/P11 );
    Lm.set(5,5,  -(P11*P22 - 4.0*P12*P12)/(3.0*P11*P11) );
    Lm.set(5,6,  -2.0*P12/P11 );
    Lm.set(5,8,   1.0 );
    
    Lm.set(6,1,  -(u1*u1*P11*P23 - 4.0*P13*u1*u1*P12 + 3.0*u1*P13*u2*P11 + 3.0*u1*P12*u3*P11 - 3.0*u2*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(6,2,   (2.0*P11*P23*u1 - 8*u1*P13*P12 + 3.0*u2*P11*P13 + 3.0*u3*P11*P12)/(3.0*P11*P11) );
    Lm.set(6,3,   (u1*P13 - u3*P11)/P11 );
    Lm.set(6,4,   (u1*P12 - u2*P11)/P11 );
    Lm.set(6,5,  -(P11*P23 - 4.0*P13*P12)/(3.0*P11*P11) );
    Lm.set(6,6,  -P13/P11 );
    Lm.set(6,7,  -P12/P11 );
    Lm.set(6,9,   1.0 );
    
    Lm.set(7,1,  -(u1*u1*P11*P33 - 4.0*u1*u1*P13*P13 + 6.0*u1*P13*u3*P11 - 3.0*u3*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(7,2,   2.0*(P11*P33*u1 - 4.0*u1*P13*P13 + 3.0*u3*P11*P13)/(3.0*P11*P11) );
    Lm.set(7,4,   (2.0*u1*P13 - 2.0*u3*P11)/P11 );
    Lm.set(7,5,  -(P11*P33 - 4.0*P13*P13)/(3.0*P11*P11) );
    Lm.set(7,7,  -2.0*P13/P11 );
    Lm.set(7,10,  1.0 );

    Lm.set(8,1,   (P11*u1*P12 - u2*P11*P11 - P12*u1*u1*c*rho + u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,2,   (-P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,3,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(8,5,  -P12/(2.0*P11*P11) );
    Lm.set(8,6,   1.0/(2.0*P11) );

    Lm.set(9,1,   (P11*P13*u1 - u3*P11*P11 - P13*u1*u1*c*rho + u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,2,  -(P11*P13 - 2.0*u1*P13*c*rho + u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,4,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(9,5,  -P13/(2.0*P11*P11) );
    Lm.set(9,7,   1.0/(2.0*P11) );
    
    Lm.set(10,1,   u1*(rho*u1*c - sq3*P11)/(6.0*c*P11*P11*rho) );
    Lm.set(10,2,   (sq3*P11 - 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(10,5,   1.0/(6.0*P11*P11) );
  
    // Project onto left eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        W.set(1,k, Lm.get(1,1)*Q.get(1,k) + Lm.get(1,2)*Q.get(2,k) + Lm.get(1,5)*Q.get(5,k) );

	W.set(2,k, Lm.get(2,1)*Q.get(1,k) + Lm.get(2,2)*Q.get(2,k) + Lm.get(2,3)*Q.get(3,k) 
	      + Lm.get(2,5)*Q.get(5,k) + Lm.get(2,6)*Q.get(6,k) );
	
	W.set(3,k, Lm.get(3,1)*Q.get(1,k) + Lm.get(3,2)*Q.get(2,k) + Lm.get(3,4)*Q.get(4,k) 
	      + Lm.get(3,5)*Q.get(5,k) + Lm.get(3,7)*Q.get(7,k) );

	W.set(4,k, Lm.get(4,1)*Q.get(1,k) + Lm.get(4,2)*Q.get(2,k) + Lm.get(4,5)*Q.get(5,k) );

	W.set(5,k, Lm.get(5,1)*Q.get(1,k) + Lm.get(5,2)*Q.get(2,k) + Lm.get(5,3)*Q.get(3,k) 
	      + Lm.get(5,5)*Q.get(5,k) + Lm.get(5,6)*Q.get(6,k) + Q.get(8,k) );

	W.set(6,k, Lm.get(6,1)*Q.get(1,k) + Lm.get(6,2)*Q.get(2,k) + Lm.get(6,3)*Q.get(3,k) 
	      + Lm.get(6,4)*Q.get(4,k) + Lm.get(6,5)*Q.get(5,k) + Lm.get(6,6)*Q.get(6,k) + 
	      + Lm.get(6,7)*Q.get(7,k) + Q.get(9,k) );

	W.set(7,k, Lm.get(7,1)*Q.get(1,k) + Lm.get(7,2)*Q.get(2,k) + Lm.get(7,4)*Q.get(4,k) 
	      + Lm.get(7,5)*Q.get(5,k) + Lm.get(7,7)*Q.get(7,k) + Q.get(10,k) );

	W.set(8,k, Lm.get(8,1)*Q.get(1,k) + Lm.get(8,2)*Q.get(2,k) + Lm.get(8,3)*Q.get(3,k) 
	      + Lm.get(8,5)*Q.get(5,k) + Lm.get(8,6)*Q.get(6,k) );

	W.set(9,k, Lm.get(9,1)*Q.get(1,k) + Lm.get(9,2)*Q.get(2,k) + Lm.get(9,4)*Q.get(4,k) 
	      + Lm.get(9,5)*Q.get(5,k) + Lm.get(9,7)*Q.get(7,k) );

	W.set(10,k, Lm.get(10,1)*Q.get(1,k) + Lm.get(10,2)*Q.get(2,k) + Lm.get(10,5)*Q.get(5,k) );
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
  
    // Set non-zero elements of matrix of left eigenvectors
    Lm.set(1,1,   u1*(sq3*P11 + rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,2,  -(sq3*P11 + 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(1,5,   1.0/(6.0*P11*P11) );
    
    Lm.set(2,1,  -(P11*u1*P12 - u2*P11*P11 + P12*u1*u1*c*rho - u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,2,   (P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(2,3,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(2,5,  -P12/(2.0*P11*P11) );
    Lm.set(2,6,   1.0/(2.0*P11) );

    Lm.set(3,1,  -(P11*P13*u1 - u3*P11*P11 + P13*u1*u1*c*rho - u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,2,   (P11*P13 + 2.0*u1*P13*c*rho - u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(3,4,  -(rho*u1*c + P11)/(2.0*c*P11*rho) );
    Lm.set(3,5,  -P13/(2.0*P11*P11) );
    Lm.set(3,7,   1.0/(2.0*P11) );
    
    Lm.set(4,1,  -(rho*u1*u1 - 3.0*P11)/(3.0*P11) );
    Lm.set(4,2,   2.0*rho*u1/(3.0*P11) );
    Lm.set(4,5,  -rho/(3.0*P11) );

    Lm.set(5,1,  -(u1*u1*P11*P22 - 4.0*u1*u1*P12*P12 + 6.0*u1*P12*u2*P11 - 3.0*u2*u2*P11*P11)/(3.0*P11*P11) );
    Lm.set(5,2,   2.0*(P11*P22*u1 - 4.0*u1*P12*P12 + 3.0*u2*P11*P12)/(3.0*P11*P11) );
    Lm.set(5,3,   (2.0*u1*P12 - 2.0*u2*P11)/P11 );
    Lm.set(5,5,  -(P11*P22 - 4.0*P12*P12)/(3.0*P11*P11) );
    Lm.set(5,6,  -2.0*P12/P11 );
    Lm.set(5,8,   1.0 );
    
    Lm.set(6,1,  -(u1*u1*P11*P23 - 4.0*P13*u1*u1*P12 + 3.0*u1*P13*u2*P11 + 3.0*u1*P12*u3*P11 - 3.0*u2*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(6,2,   (2.0*P11*P23*u1 - 8*u1*P13*P12 + 3.0*u2*P11*P13 + 3.0*u3*P11*P12)/(3.0*P11*P11) );
    Lm.set(6,3,   (u1*P13 - u3*P11)/P11 );
    Lm.set(6,4,   (u1*P12 - u2*P11)/P11 );
    Lm.set(6,5,  -(P11*P23 - 4.0*P13*P12)/(3.0*P11*P11) );
    Lm.set(6,6,  -P13/P11 );
    Lm.set(6,7,  -P12/P11 );
    Lm.set(6,9,   1.0 );
    
    Lm.set(7,1,  -(u1*u1*P11*P33 - 4.0*u1*u1*P13*P13 + 6.0*u1*P13*u3*P11 - 3.0*u3*u3*P11*P11)/(3.0*P11*P11) );
    Lm.set(7,2,   2.0*(P11*P33*u1 - 4.0*u1*P13*P13 + 3.0*u3*P11*P13)/(3.0*P11*P11) );
    Lm.set(7,4,   (2.0*u1*P13 - 2.0*u3*P11)/P11 );
    Lm.set(7,5,  -(P11*P33 - 4.0*P13*P13)/(3.0*P11*P11) );
    Lm.set(7,7,  -2.0*P13/P11 );
    Lm.set(7,10,  1.0 );

    Lm.set(8,1,   (P11*u1*P12 - u2*P11*P11 - P12*u1*u1*c*rho + u1*u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,2,   (-P11*P12 + 2.0*u1*P12*c*rho - u2*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(8,3,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(8,5,  -P12/(2.0*P11*P11) );
    Lm.set(8,6,   1.0/(2.0*P11) );

    Lm.set(9,1,   (P11*P13*u1 - u3*P11*P11 - P13*u1*u1*c*rho + u1*u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,2,  -(P11*P13 - 2.0*u1*P13*c*rho + u3*c*P11*rho)/(2.0*c*P11*P11*rho) );
    Lm.set(9,4,   (P11 - rho*u1*c)/(2.0*c*P11*rho) );
    Lm.set(9,5,  -P13/(2.0*P11*P11) );
    Lm.set(9,7,   1.0/(2.0*P11) );
    
    Lm.set(10,1,   u1*(rho*u1*c - sq3*P11)/(6.0*c*P11*P11*rho) );
    Lm.set(10,2,   (sq3*P11 - 2.0*rho*u1*c)/(6.0*c*P11*P11*rho) );
    Lm.set(10,5,   1.0/(6.0*P11*P11) );    
  
    // Project onto left eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        W.set(11,k, Lm.get(1,1)*Q.get(11,k) + Lm.get(1,2)*Q.get(12,k) + Lm.get(1,5)*Q.get(15,k) );

	W.set(12,k, Lm.get(2,1)*Q.get(11,k) + Lm.get(2,2)*Q.get(12,k) + Lm.get(2,3)*Q.get(13,k) 
	      + Lm.get(2,5)*Q.get(15,k) + Lm.get(2,6)*Q.get(16,k) );
	
	W.set(13,k, Lm.get(3,1)*Q.get(11,k) + Lm.get(3,2)*Q.get(12,k) + Lm.get(3,4)*Q.get(14,k) 
	      + Lm.get(3,5)*Q.get(15,k) + Lm.get(3,7)*Q.get(17,k) );

	W.set(14,k, Lm.get(4,1)*Q.get(11,k) + Lm.get(4,2)*Q.get(12,k) + Lm.get(4,5)*Q.get(15,k) );

	W.set(15,k, Lm.get(5,1)*Q.get(11,k) + Lm.get(5,2)*Q.get(12,k) + Lm.get(5,3)*Q.get(13,k) 
	      + Lm.get(5,5)*Q.get(15,k) + Lm.get(5,6)*Q.get(16,k) + Q.get(18,k) );

	W.set(16,k, Lm.get(6,1)*Q.get(11,k) + Lm.get(6,2)*Q.get(12,k) + Lm.get(6,3)*Q.get(13,k) 
	      + Lm.get(6,4)*Q.get(14,k) + Lm.get(6,5)*Q.get(15,k) + Lm.get(6,6)*Q.get(16,k) + 
	      + Lm.get(6,7)*Q.get(17,k) + Q.get(19,k) );

	W.set(17,k, Lm.get(7,1)*Q.get(11,k) + Lm.get(7,2)*Q.get(12,k) + Lm.get(7,4)*Q.get(14,k) 
	      + Lm.get(7,5)*Q.get(15,k) + Lm.get(7,7)*Q.get(17,k) + Q.get(20,k) );

	W.set(18,k, Lm.get(8,1)*Q.get(11,k) + Lm.get(8,2)*Q.get(12,k) + Lm.get(8,3)*Q.get(13,k) 
	      + Lm.get(8,5)*Q.get(15,k) + Lm.get(8,6)*Q.get(16,k) );

	W.set(19,k, Lm.get(9,1)*Q.get(11,k) + Lm.get(9,2)*Q.get(12,k) + Lm.get(9,4)*Q.get(14,k) 
	      + Lm.get(9,5)*Q.get(15,k) + Lm.get(9,7)*Q.get(17,k) );

	W.set(20,k, Lm.get(10,1)*Q.get(11,k) + Lm.get(10,2)*Q.get(12,k) + Lm.get(10,5)*Q.get(15,k) );
    }
    

    // -----------------------------
    // For the ELECTROMAGNETIC FIELD
    // -----------------------------
    
    // Project onto left eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        W.set(21,k, (cs_light*Q.get(22,k) + Q.get(26,k))/(2.0*cs_light) );

	W.set(22,k, (cs_light*Q.get(23,k) - Q.get(25,k))/(2.0*cs_light) );

	W.set(23,k, Q.get(21,k) );
	
	W.set(24,k, Q.get(24,k) );

	W.set(25,k, (cs_light*Q.get(22,k) - Q.get(26,k))/(2.0*cs_light) );

	W.set(26,k, (cs_light*Q.get(23,k) + Q.get(25,k))/(2.0*cs_light) );
    }
}
