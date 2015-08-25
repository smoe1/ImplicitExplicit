/* poly/solve_cubic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */

#include <math.h>

// generic swap that uses gcc __typeof__ extension
//#define Swap(X,Y)  do{ __typeof__ (X) _T = X; X = Y; Y = _T; }while(0)
#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

/* I rewrote this to eliminate division (in case not using fast math). -eaj
   */
int 
/*
gsl_poly_solve_cubic (double a, double b, double c, 
                      double *x0, double *x1, double *x2)
*/
gsl_poly_solve_cubic (const double* p, double *roots)
{
  const double& c = p[0];
  const double& b = p[1];
  const double& a = p[2];
  const double a_3 = a*(1./3.);
  double* x0 = &roots[0];
  double* x1 = &roots[1];
  double* x2 = &roots[2];

  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double CR2 = 729 * r * r; // 729 = 9*9*9
  double CQ3 = 2916 * q * q * q; // 2916 = 54*54

  double Q = q *(1./9.);
  double R = r *(1./54.);

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  if (R == 0 && Q == 0)
    {
      *x0 = - a_3 ;
      *x1 = - a_3 ;
      *x2 = - a_3 ;
      return 3 ;
    }
  else if (CR2 == CQ3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          *x0 = -2 * sqrtQ  - a_3;
          *x1 = sqrtQ - a_3;
          *x2 = sqrtQ - a_3;
        }
      else
        {
          *x0 = - sqrtQ  - a_3;
          *x1 = - sqrtQ - a_3;
          *x2 = 2 * sqrtQ - a_3;
        }
      return 3 ;
    }
  else if (R2 < Q3)
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double ratio = sgnR * sqrt (R2 / Q3);
      double theta = acos (ratio); // very expensive
      double norm = -2 * sqrt (Q);
      /* cos is expensive.  could use angle addition formulas
         (which still involves evaluation of sin and cos) -eaj */
      *x0 = norm * cos (theta * (1./3)) - a_3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) *(1./3.)) - a_3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) *(1./3.)) - a_3;
      
      /* Sort *x0, *x1, *x2 into increasing order */

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;
      
      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;
          
          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }
      
      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      /* I changed pow(*,1./3.) to cbrt. -eaj */
      double A = -sgnR * cbrt(fabs (R) + sqrt (R2 - Q3));
      double B = Q / A ;
      double ma_3 = -a_3;
      double A_plus_B = A + B;
      /* eaj: could modify this to return the real part of
       * the complex conjugate pair for the other two roots.
       * I think that what I am doing is a close approximation to that. */
      /* eaj: set the other roots to be the extremum most distant
         from the root; else the point of inflection.
         This handles the case of missing real roots due
         to finite precision. */
      /* sorting the "roots" eliminates backward compatibility 
       * with the original GSL implementation */
      #define SORT_ROOTS 1
      double true_root = A_plus_B + ma_3;
      double disc = ma_3*ma_3 - b*(1./3.);
      // use the extremum farthest from the root
      //
      double sqd = 0.; // if no extremum then use inflection point
      if(disc >=0) sqd = sqrt(disc);
      if(A_plus_B > 0)
      {
        #ifdef SORT_ROOTS
        *x0 = *x1 = ma_3 - sqd;
        *x2 = true_root;
        #else
        *x0 = true_root;
        *x1 = *x2 = ma_3 - sqd;
        #endif
      }
      else
      {
        *x0 = true_root;
        *x1 = *x2 = ma_3 + sqd;
      }
      return 1;
    }
}
