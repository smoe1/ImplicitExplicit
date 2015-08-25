#include "dogdefs.h"
#include "dog_math.h"

double GetAlpha(const double rho,
                const double p,
                const double q,
                const double r)
{  
  // check inputs
  const double p2 = p*p;
  const double p3 = p*p2;
  const double q2 = q*q;
  const double rmin = (p2/rho + q2/p) - 1.0e-8;
  const double rmax = 3.0*p2/rho + 1.0e-8;
  /*
  if (rho<=0.0 || p<=0.0 || r<rmin || ((fabs(q)<=1.0e-14) && r>rmax) )
    {
      printf("  rho = %22.16e \n",rho);
      printf("    p = %22.16e \n",p);
      printf("    q = %22.16e \n",q);
      printf("    r = %22.16e \n",r);
      printf("\n");
      printf(" rmin = %22.16e \n",rmin);
      printf(" rmax = %22.16e \n",rmax);
      printf("\n");
      exit(1);
    }
  */
  // compute alpha
  if (fabs(q) <= 1.0e-12)
    {
      double tmp = 3.0*p2-r*rho;
      if (tmp > 0.0)
	{
	  return sqrt(tmp/(2.0*p2));          
	}
      else
	{
	  return 0.0;
     	}
    }
  else
    { 
      double Q = 0.5*onethird*(3.0-rho*r/p2);
      double R = rho*q2/(4.0*p3);
      double R2 = R*R;
      double Q3 = pow(Q,3);
      double D = R2-Q3;
      double SD = sqrt(D);
      double SQ = sqrt(Q);
      
      if ( D >= 0.0 )
        {
          double T;
          double Tstar;
          Tstar = R - SD;
          if ( Tstar < 0.0 )
            {
              T = -pow(fabs(Tstar),onethird);
            }
          else 
            {
              T = pow(fabs(Tstar),onethird);
            }
          return pow(R+SD,onethird)+T;       
        }
      else
        {
          double theta = onethird * acos(R/sqrt(Q3));
          double a1 = 2.0*SQ*cos(theta);
          if ( a1>=0.0 && a1<=1.0 )
            {
              return a1;                         
            }
          else
            {
              double a2 = 2.0*SQ*cos(theta+2.0*onethird*pi);
              if ( a2>=0.0 && a2<=1.0 )
                {
                  return a2;                                
                }
              else
                {
                  double a3 = 2.0*SQ*cos(theta+4.0*onethird*pi);
                  return a3;                 
                }
            }
        }
    }

}
