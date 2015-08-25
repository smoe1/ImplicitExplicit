#include "dogdefs.h"
#include "dog_math.h"

void GetWeights(const double rho, const double Qt123, double mut[][3], double* w)
{
   double A[8][8];
   double b[8] = {rho,0.0,0.0,0.0,0.0,0.0,0.0,Qt123};
   double U[8][8];
   double L[8][8];
   double y[8];
   double p = 0.0, temp;
   int    PI[8];
   int    k;

   for (int i=0; i<8; i++)
     {
        A[0][i] = 1.0;
        A[1][i] = mut[i][0];
        A[2][i] = mut[i][1];
        A[3][i] = mut[i][2];
        A[4][i] = mut[i][0]*mut[i][1];
        A[5][i] = mut[i][0]*mut[i][2];
        A[6][i] = mut[i][1]*mut[i][2];
        A[7][i] = mut[i][0]*mut[i][1]*mut[i][2];

        L[i][i] = 1.0;
        for (int j=i+1; j<8; j++)
          {
             L[i][j] = 0.0;
          }
    
        for (int j=0;j<i;j++)
          {
             U[i][j] = 0.0;
          }

        PI[i] = i;
     }

   for (int i=0; i<8; i++)  
     {
        p = 0.0;
        k = i;
        for(int j=i; j<8; j++)
           {
              if (abs(A[j][i]) > p)
                 {
                    p = abs(A[j][i]);
                    k = j;
                 } 
           }
     
         if (p==0.0)
           {
              printf("ERROR: Singular\n");
              printf("\n");
              exit(1);
           }

         temp  = PI[i];
         PI[i] = PI[k];
         PI[k] = temp;
  
         for (int j=0;j<8;j++)
           {
              temp = A[i][j];
              A[i][j] = A[k][j];
              A[k][j] = temp;
           }


         for (int j=i+1;j<8;j++)
           {
              A[j][i] = A[j][i]/A[i][i];
           }

         for (int j=i+1;j<8;j++)
           {
              for (int s=i+1;s<8;s++)
                 {
                    A[j][s] = A[j][s] - A[j][i]*A[i][s];
                 }
           }
      }

    for (int i=0;i<8;i++)
      {
         for (int j=0;j<i-1;j++)
           {
              L[i][j] = A[i][j];
           }
         for (int j=i;j<8;j++)
           {
              U[i][j] = A[i][j];
           }
      }

    for (int i=0;i<8;i++)
      {
         y[i] = b[PI[i]];
         
         for (int j=0;j<i-1;j++)
           { 
             y[i] -= L[i][j]*y[j];
           }
      }

    for (int i=7;i>=0;i--)
      {
         w[i] = y[i];
         
         for (int j=i+1;j<8;j++)
           {
             w[i] -= U[i][j]*w[j];
           }
         w[i]/=U[i][i];
      }
}  