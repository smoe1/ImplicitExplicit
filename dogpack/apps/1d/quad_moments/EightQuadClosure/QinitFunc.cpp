#include "dogdefs.h"
#include <cmath>
#include "QuadMomentParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts = xpts.getsize();

    // density, velocity, pressure and (heat?) tensor
    double rho,u1,u2,u3,P11,P22,P33,P12,P13,P23,
           Q111,Q112,Q113,Q122,Q123,Q133,Q222,Q223,Q233,Q333;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        /*

        // uniform distribution       

<<<<<<< HEAD
  // uniform distribution       
    
      rho = 2.0;
      u1  = 2.0;
      u2  = 2.0;
      u3  = 2.0;
      P11 = 2.0;
      P22 = 2.0;
      P33 = 2.0;
      P12 = 0.0;
      P13 = 0.0;
      P23 = 0.0;
      Q111= 2.0;
      //Q112= 0.6;
      //Q113= 0.6;
      //Q122= 0.6;
      Q123= 0.6;
      //Q133= 0.6;
      Q222= 2.0;
      //Q223= 0.6;
      //Q233= 0.6;
      Q333= 2.0;
 */
   // shock-tube
=======
        rho = 2.0;
        u1  = 2.0;
        u2  = 2.0;
        u3  = 2.0;
        P11 = 2.0;
        P22 = 2.0;
        P33 = 2.0;
        P12 = 0.0;
        P13 = 0.0;
        P23 = 0.0;
        Q111= 2.0;
        Q112= 0.6;
        Q113= 0.6;
        Q122= 0.6;
        Q123= 0.6;
        Q133= 0.6;
        Q222= 2.0;
        Q223= 0.6;
        Q233= 0.6;
        Q333= 2.0;
         */
        // shock-tube
>>>>>>> d1b88df4d167b46dbe8b35ed7cbda2424b0ac117

        if (x<0.5e0)
        {
            rho = 1.0;
            u1  = 0.0;
            u2  = 0.0;
            u3  = 0.0;
            P11 = 1.0;
            P22 = 2.0;
            P33 = 3.0;
            P12 = 0.25;
            P13 = 0.0;
            P23 = 0.0;
            Q111= 0.0;
            //Q112= 0.0;
            //Q113= 0.0;
            //Q122= 0.0;
            Q123= 0.0;
            //Q133= 0.0;
            Q222= 0.0;
            //Q223= 0.0;
            //Q233= 0.0;
            Q333= 0.0;
        }
        else
        {
            rho = 2.0;
            u1  = 0.0;
            u2  = 0.0;
            u3  = 0.0;
            P11 = 2.0;
            P22 = 2.0;
            P33 = 2.0;
            P12 = 0.0;
            P13 = 0.0;
            P23 = 0.0;
            Q111= 0.0;
            //Q112= 0.0;
            //Q113= 0.0;
            //Q122= 0.0;
            Q123= 0.0;
<<<<<<< HEAD
            //Q133= 0.0;
            Q222= 0.0;
            //Q223= 0.0;
            //Q233= 0.0;
            Q333= 0.0;
         }
  
                                                              
      // Check input
      const double detP = P11*P22*P33-P11*pow(P23,2)-P22*pow(P13,2)-P33*pow(P12,2)+2.0*P12*P13*P23;
      if (rho<=0.0 || P11<=0.0 || P22<=0.0 || P33<=0.0 || P11*P22-pow(P12,2)<=0.0 || detP<=0.0)
=======
            Q133= 0.0;
            Q222= 2.0;
            Q223= 0.0;
            Q233= 0.0;
            Q333= 2.0;
        }


        // Check input
        const double detP = P11*P22*P33-P11*pow(P23,2)-P22*pow(P13,2)-P33*pow(P12,2)+2.0*P12*P13*P23;
        if (rho<=0.0 || P11<=0.0 || P22<=0.0 || P33<=0.0 || P11*P22-pow(P12,2)<=0.0 || detP<=0.0)
>>>>>>> d1b88df4d167b46dbe8b35ed7cbda2424b0ac117
        {
            printf("              rho = %22.16e \n",rho);
            printf("              P11 = %22.16e \n",P11);
            printf("              P22 = %22.16e \n",P22);
            printf("              P33 = %22.16e \n",P33);
            printf("  P11*P22-P12*P12 = %22.16e \n",P11*P22-pow(P12,2));
            printf("             detP = %22.16e \n",detP);
            printf("\n");
            exit(1);
        }

<<<<<<< HEAD
      const double M000 = rho;
      const double M100 = rho*u1;
      const double M010 = rho*u2;
      const double M001 = rho*u3;
      const double M200 = P11+rho*u1*u1;
      const double M020 = P22+rho*u2*u2;
      const double M002 = P33+rho*u3*u3;
      const double M110 = P12+rho*u1*u2;
      const double M101 = P13+rho*u1*u3;
      const double M011 = P23+rho*u2*u3;
      const double M300 = Q111+3.0*u1*P11+rho*u1*u1*u1;
      //const double M210 = Q112+2.0*u1*P12+u2*P11+rho*u1*u1*u2;
      //const double M201 = Q113+2.0*u1*P13+u3*P11+rho*u1*u1*u3;
      //const double M120 = Q122+2.0*u2*P12+u1*P22+rho*u1*u2*u2;
      const double M111 = Q123+u1*P23+u2*P13+u3*P12+rho*u1*u2*u3;
      //const double M102 = Q133+2.0*u3*P13+u1*P33+rho*u1*u3*u3;
      const double M030 = Q222+3.0*u2*P22+rho*u2*u2*u2;
      //const double M021 = Q223+2.0*u2*P23+u3*P22+rho*u2*u2*u3;
      const double M003 = Q333+3.0*u3*P33+rho*u3*u3*u3;
      
      /* 
      qvals.set(i, 1, M000  );  // density
      qvals.set(i, 2, M100  );  // momentum
      qvals.set(i, 3, M010  );  
      qvals.set(i, 4, M001  );  
      qvals.set(i, 5, M200  );  // energy
      qvals.set(i, 6, M110  );
      qvals.set(i, 7, M101  );
      qvals.set(i, 8, M020  );
      qvals.set(i, 9, M011  );
      qvals.set(i,10, M002  );
      qvals.set(i,11, M300  );  // heat flow
      qvals.set(i,12, M210  );
      qvals.set(i,13, M201  );
      qvals.set(i,14, M120  );
      qvals.set(i,15, M111  );
      qvals.set(i,16, M102  );
      qvals.set(i,17, M030  );
      qvals.set(i,18, M021  );
      qvals.set(i,19, M021  );
      qvals.set(i,20, M003  );  */

      qvals.set(i, 1, M000  );  // density
      qvals.set(i, 2, M100  );  // momentum
      qvals.set(i, 3, M010  );  
      qvals.set(i, 4, M001  );  
      qvals.set(i, 5, M200  );  // energy
      qvals.set(i, 6, M110  );
      qvals.set(i, 7, M101  );
      qvals.set(i, 8, M020  );
      qvals.set(i, 9, M011  );
      qvals.set(i,10, M002  );  // heat flow
      qvals.set(i,11, M300  );
      qvals.set(i,12, M030  );
      qvals.set(i,13, M003  );
      qvals.set(i,14, M111  );
      
=======
        const double M000 = rho;
        const double M100 = rho*u1;
        const double M010 = rho*u2;
        const double M001 = rho*u3;
        const double M200 = P11+rho*u1*u1;
        const double M020 = P22+rho*u2*u2;
        const double M002 = P33+rho*u3*u3;
        const double M110 = P12+rho*u1*u2;
        const double M101 = P13+rho*u1*u3;
        const double M011 = P23+rho*u2*u3;
        const double M300 = Q111+3.0*u1*P11+rho*u1*u1*u1;
        const double M210 = Q112+2.0*u1*P12+u2*P11+rho*u1*u1*u2;
        const double M201 = Q113+2.0*u1*P13+u3*P11+rho*u1*u1*u3;
        const double M120 = Q122+2.0*u2*P12+u1*P22+rho*u1*u2*u2;
        const double M111 = Q123+u1*P23+u2*P13+u3*P12+rho*u1*u2*u3;
        const double M102 = Q133+2.0*u3*P13+u1*P33+rho*u1*u3*u3;
        const double M030 = Q222+3.0*u2*P22+rho*u2*u2*u2;
        const double M021 = Q223+2.0*u2*P23+u3*P22+rho*u2*u2*u3;
        const double M003 = Q333+3.0*u3*P33+rho*u3*u3*u3;

        qvals.set(i, 1, M000  );  // density
        qvals.set(i, 2, M100  );  // momentum
        qvals.set(i, 3, M010  );  
        qvals.set(i, 4, M001  );  
        qvals.set(i, 5, M200  );  // energy
        qvals.set(i, 6, M110  );
        qvals.set(i, 7, M101  );
        qvals.set(i, 8, M020  );
        qvals.set(i, 9, M011  );
        qvals.set(i,10, M002  );
        qvals.set(i,11, M300  );  // heat flow
        qvals.set(i,12, M210  );
        qvals.set(i,13, M201  );
        qvals.set(i,14, M120  );
        qvals.set(i,15, M111  );
        qvals.set(i,16, M102  );
        qvals.set(i,17, M030  );
        qvals.set(i,18, M021  );
        qvals.set(i,19, M021  );
        qvals.set(i,20, M003  );
>>>>>>> d1b88df4d167b46dbe8b35ed7cbda2424b0ac117
    }

}
