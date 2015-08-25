#include "dogdefs.h"
#include "dog_math.h"

void LinTrans(const double* Q, const double* L, double* Qt)

{
   const double L11  = L[0];
   const double L22  = L[1];
   const double L33  = L[2];
   const double L21  = L[3];
   const double L31  = L[4];
   const double L32  = L[5];
   /* 
   const double Q111 = Q[0];
   const double Q112 = Q[1];
   const double Q113 = Q[2];
   const double Q122 = Q[3];
   const double Q123 = Q[4];
   const double Q133 = Q[5];
   const double Q222 = Q[6];
   const double Q223 = Q[7];
   const double Q233 = Q[8];
   const double Q333 = Q[9];  */
  
   const double Q111 = Q[0];
   const double Q222 = Q[1];
   const double Q333 = Q[2];
   const double Q123 = Q[3];

 

   //Find inverse of L
   const double T11 = 1/L11;
   const double T22 = 1/L22;
   const double T33 = 1/L33;
   const double T21 = -L21/L11/L22;
   const double T31 = (L21*L32-L22*L31)/(L11*L22*L33);
   const double T32 = -L32/L22/L33;


   //Convert Q to Qt
/*
   Qt[0] = pow(T11,3)*Q111;
   Qt[1] = pow(T11,2)*T21*Q111 + pow(T11,2)*T22*Q112;
   Qt[2] = pow(T11,2)*T31*Q111 + pow(T11,2)*T32*Q112 + pow(T11,2)*T33*Q113;
   Qt[3] = T11*pow(T21,2)*Q111 + 2.0*T11*T21*T22*Q112 + T11*pow(T22,2)*Q122;
   Qt[4] = T11*T21*T31*Q111 + T11*(T21*T32+T22*T31)*Q112 
                         + T11*T21*T33*Q113 + T11*T22*T32*Q122 + T11*T22*T33*Q123;
   Qt[5] = T11*pow(T31,2)*Q111 + 2.0*T11*T31*T32*Q112 + 2.0*T11*T31*T33*Q113
                         + T11*pow(T32,2)*Q122 + 2.0*T11*T32*T33*Q123 + T11*pow(T33,2)*Q133;
   Qt[6] = pow(T21,3)*Q111 + 3.0*T21*T21*T22*Q112 + 3.0*T21*T22*T22*Q122 + pow(T22,3)*Q222;
   Qt[7] = pow(T21,2)*T31*Q111 + T21*(T21*T32+2.0*T22*T31)*Q112 + pow(T21,2)*T33*Q113
                         + T22*(2.0*T21*T32+T22*T31)*Q122 + 2.0*T21*T22*T33*Q123 
                         + pow(T22,2)*T32*Q222 + pow(T22,2)*T33*Q223;
   Qt[8] = T21*pow(T31,2)*Q111 + T31*(2.0*T21*T32+T22*T31)*Q112 + 2.0*T21*T31*T33*Q113
                         + 2.0*T21*T32*T33*Q123 + T21*pow(T33,2)*Q133 + T32*(T21*T32+2.0*T22*T31)*Q122
                         + 2.0*T22*T31*T33*Q123 + T22*pow(T32,2)*Q222 + 2.0*T22*T32*T33*Q223 + T22*pow(T33,2)*Q233;
   Qt[9] = pow(T31,3)*Q111 + 3.0*pow(T31,2)*T32*Q112 + 3.0*pow(T31,2)*T33*Q113 + 3.0*T31*pow(T32,2)*Q122
                         + 6.0*T31*T32*T33*Q123 + 3.0*T31*pow(T33,2)*Q133 + T32*pow(T32,2)*Q222 + 3.0*pow(T32,2)*T33*Q223
                         + 3.0*T32*pow(T33,2)*Q233 + pow(T33,3)*Q333;
*/

   Qt[0] = pow(T11,3)*Q111;  //Qt111
   Qt[1] = pow(T21,3)*Q111+pow(T22,3)*Q222;   //Qt222
   Qt[2] = pow(T32,3)*Q222+(pow(T31,3)+3.0*T21*pow(T31,2)*T32/T22-9.0*pow(T21,2)*T31*pow(T32,2)/pow(T22,2)+6.0*pow(T21,3)*pow(T32,3)/pow(T22,3))*Q111+pow(T33,3)*Q333+(6.0*T21*pow(T32,2)*T33/T22-6.0*T31*T32*T33)*Q123;   //Qt333
   Qt[3] = (T11*pow(T21,2)*T32/T22-T11*T21*T31)*Q111+T11*T22*T33*Q123;   //Qt123

}