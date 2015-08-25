#include "dogdefs.h"
#include "math.h"
#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "DogParams.h"
#include <stdlib.h>
#include <string>
#include "QuadMomentParams.h"

// Function that is called before each time step
void BeforeStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   mbc  = q.getmbc();
  const int   maux = aux.getsize(2);
  
  void L2Project(int mopt, int istart, int iend,
                 const dTensor2& node,
                 const dTensorBC3& qin, 
                 const dTensorBC3& auxin,  
                 dTensorBC3& Fout,
                 void (*Func)(const dTensor1&, const dTensor2&, 
                              const dTensor2&, dTensor2&));
  void Rootfinding(const dTensor1& xpts, 
                   const dTensor2& Q, 
                   const dTensor2& Aux, 
                   dTensor2& auxvals);

  const string closure = quadMomentParams.closure;

  int cl=0;
  if (closure == "M4")
    {  cl = 1;  }
  else if (closure == "M6")
    {  cl = 2;  }
  else if (closure == "M8")
    {  cl = 3;  }

  
  switch(cl){
    case 1:
       break;
    case 2:
       L2Project(0,1-mbc,melems+mbc,node,q,aux,aux,&Rootfinding);
       break;
    case 3:
       L2Project(0,1-mbc,melems+mbc,node,q,aux,aux,&Rootfinding);
       break;
    }
       
}
   

void BeforeStep(double dt, const dTensor2& node, dTensorBC3& aux, dTensorBC3& q, 
                    void* data)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   maux = aux.getsize(2);  
}


void Rootfinding(const dTensor1& xpts, 
                   const dTensor2& Q, 
                   const dTensor2& Aux, 
                   dTensor2& auxvals)
{
  const int numpts=xpts.getsize();
  const string closure = quadMomentParams.closure;
  double* x;
  double* y;
  double Laguerres(double a[],int n,double z);
  double BirgeVieta(double a[], int n, double z);

  int cl = 0;
  if (closure == "M4")
    {  cl = 1;  }
  else if (closure == "M6")
    {  cl = 2;  }
  else if (closure == "M8")
    {  cl = 3;  }
  
  switch(cl)
  {
    case 1:
   //////////////////////////////// M4 closure /////////////////////////////

    break;

    case 2:
   //////////////////////////////// M6 closure //////////////////////////////
    {
     x = new double[3];  
     y = new double[3];

    for (int i=1; i<=numpts; i++)
      {  
        const double rho    = Q.get(i,1);
        const double m1     = Q.get(i,2);
        const double u      = Q.get(i,2)/rho;
	const double m2     = Q.get(i,3);
	const double press  = Q.get(i,3)-rho*pow(u,2);
	const double q      = Q.get(i,4)-rho*pow(u,3)-3.0*press*u;  
	const double r      = Q.get(i,5)-4.0*q*u-6.0*press*pow(u,2)-rho*pow(u,4);
	const double s      = Q.get(i,6)-5.0*r*u-10.0*q*pow(u,2)-10.0*press*pow(u,3)-rho*pow(u,5); 
        printf("rho = %f\n",rho);
        printf("u = %f\n",u);
        printf("press = %f\n",press);
        printf("q = %f\n",q);
        printf("r = %f\n",r);
        printf("s = %f\n",s);
	  
	const double a = rho*pow(q,2)+pow(press,3)-rho*press*r;
	const double b = rho*press*s-q*pow(press,2)-rho*q*r;
	const double c = rho*pow(r,2)+press*pow(q,2)-rho*s*q-r*pow(press,2);
	const double d = 2.0e0*press*q*r-pow(q,3)-s*pow(press,2);
        
        double p1 = -1.0/3.0*pow(b/a,2)+c/a;
        double p2 = 2.0/27.0*pow(b/a,3)-b*c/3.0/pow(a,2)+d/a;

        const int MaxIters = 1000;
        const double TOL = 1.0e-12;
        const double disc = pow(p1/3.0,3)+pow(p2/2.0,2);

        x[0] = Aux.get(i,1)-u;
        x[1] = Aux.get(i,2)-u;
        x[2] = Aux.get(i,3)-u;

        if(disc<=1.0e-12)
          {
            int mstop = 0;
            int NumIters = 0;
            
            while (mstop==0)
              {
                 double f[3] = {x[0]+x[1]+x[2]+b/a, 
			 x[0]*x[1]+x[0]*x[2]+x[1]*x[2]-c/a, 
			 x[0]*x[1]*x[2]+d/a};
                 double invJ[3][3] = {pow(x[0],2)/(x[0]-x[1])/(x[0]-x[2]),
			       -x[0]/(x[0]-x[1])/(x[0]-x[2]),
			       1.0/(x[0]-x[1])/(x[0]-x[2]),
			       -pow(x[1],2)/(x[0]-x[1])/(x[1]-x[2]),
			       x[1]/(x[0]-x[1])/(x[1]-x[2]),
			       -1.0/(x[0]-x[1])/(x[1]-x[2]),
			       pow(x[2],2)/(x[0]-x[2])/(x[1]-x[2]),
			       -x[2]/(x[0]-x[2])/(x[1]-x[2]),
			       1.0/(x[0]-x[2])/(x[1]-x[2])};
	         double xnew[3];
                 for (int i=0; i<3; i++)
	           {
	              xnew[i] = x[i] - (invJ[i][0]*f[0]+invJ[i][1]*f[1]+invJ[i][2]*f[2]);
	           }
	  
	         double error = Max(Max(fabs(xnew[0]-x[0]),fabs(xnew[1]-x[1])),fabs(xnew[2]-x[2]));
	         x[0] = xnew[0];
	         x[1] = xnew[1];
	         x[2] = xnew[2];  
	         NumIters = NumIters + 1;
                 
                 if (error <= TOL || NumIters == MaxIters)
	           { 
	              mstop = 1; 
	              printf(" NumIters = %i\n",NumIters);
	           }
	  
	      }
           }
         else
           {
                 printf(" ERROR in RootFinding: Not all roots are real\n");
                 printf("       a = %e\n",a);
                 printf("       b = %e\n",b);
                 printf("       c = %e\n",c);
                 printf("       d = %e\n",d);
                 printf("       p1 = %e\n",p1);
                 printf("       p2 = %e\n",p2);
                 printf("    disc = %e\n",disc);
                 exit(1);
           }
      y[0] = (m2-m1*(x[1]+u)-m1*(x[2]+u)+rho*(x[1]+u)*(x[2]+u))/((x[0]-x[1])*(x[0]-x[2]));
      y[1] = (-m2+m1*(x[0]+u)+m1*(x[2]+u)-rho*(x[0]+u)*(x[2]+u))/((x[0]-x[1])*(x[1]-x[2]));
      y[2] = (m2-m1*(x[0]+u)-m1*(x[1]+u)+rho*(x[0]+u)*(x[1]+u))/((x[0]-x[2])*(x[1]-x[2]));

      auxvals.set(i,1,x[0]+u);
      auxvals.set(i,2,x[1]+u);
      auxvals.set(i,3,x[2]+u);
      auxvals.set(i,4,y[0]);
      auxvals.set(i,5,y[1]);
      auxvals.set(i,6,y[2]);

     }
    break;
   }

    case 3:

 /////////////////////////////////// M8 closure ////////////////////////////////////
    {
     x = new double[4];
     y = new double[4];
     double a[5];
     double b[4];
     double c[3];

    for (int i=1;i<=numpts;i++)
      {
        const double rho    = Q.get(i,1);
        const double m1     = Q.get(i,2);
        const double u      = Q.get(i,2)/rho;
	const double m2     = Q.get(i,3);
	const double press  = Q.get(i,3)-rho*pow(u,2);
        const double m3     = Q.get(i,4);
	const double q      = Q.get(i,4)-rho*pow(u,3)-3.0*press*u; 
        const double m4     = Q.get(i,5); 
	const double r      = Q.get(i,5)-4.0*q*u-6.0*press*pow(u,2)-rho*pow(u,4);
	const double s      = Q.get(i,6)-5.0*r*u-10.0*q*pow(u,2)
                              -10.0*press*pow(u,3)-rho*pow(u,5); 
        const double delta  = Q.get(i,7)-6.0*s*u-15.0*r*pow(u,2)
                              -20.0*q*pow(u,3)-15.0*press*pow(u,4)-rho*pow(u,6);
        const double tau    = Q.get(i,8)-7.0*delta*u-21.0*s*pow(u,2)-35.0*r*pow(u,3)
                              -35.0*q*pow(u,4)-21.0*press*pow(u,5)-rho*pow(u,7);
/*
        printf("rho = %f\n",rho);
        printf("u = %f\n",u);
        printf("press = %f\n",press);
        printf("q = %f\n",q);
        printf("r = %f\n",r);
        printf("s = %f\n",s);
        printf("delta = %f\n",delta);
        printf("tau = %f\n",tau);*/


        a[4] = delta*pow(press,3)-pow(r*press,2)-2.0*s*q*pow(press,2)
                         +3.0*r*press*pow(q,2)+press*pow(s,2)*rho-press*delta*rho*r
                         -pow(q,4)+delta*rho*pow(q,2)+rho*pow(r,3)-2.0*rho*r*s*q;
        a[3] = -rho*r*delta*q+press*s*rho*delta-press*tau*rho*r
                         -pow(press,2)*delta*q-r*pow(q,3)+pow(press,3)*tau
                         -2.0*pow(press,2)*s*r+s*pow(r,2)*rho-pow(s,2)*rho*q
                         +tau*rho*pow(q,2)+2.0*pow(r,2)*q*press+press*s*pow(q,2);
        a[2] = -pow(r*q,2)+s*pow(q,3)+press*pow(r,3)-q*rho*s*delta+r*rho*pow(s,2)
                         -pow(r,2)*rho*delta+press*rho*pow(delta,2)-press*pow(q,2)*delta
                         -pow(press,2)*r*delta+q*r*rho*tau-press*rho*s*tau
                         +pow(press,2)*q*tau;
        a[1] = -q*pow(r,3)+rho*pow(s,3)+pow(press,2)*s*delta
                         -pow(q,3)*delta-2.0*press*q*pow(s,2)
                         +press*pow(q,2)*tau-2.0*r*rho*s*delta-q*rho*s*tau
                         +press*s*pow(r,2)+pow(r,2)*rho*tau+q*rho*pow(q,2)
                         -pow(press,2)*r*tau;
        a[0] = -pow(s*q,2)+pow(q,3)*tau-pow(press*delta,2)-pow(r,4)
                         -2.0*pow(q,2)*r*delta-2.0*press*r*pow(s,2)
                         +2.0*press*pow(r,2)*delta+pow(press,2)*s*tau
                         +3.0*q*pow(r,2)*s+2.0*press*q*s*delta-2.0*press*q*r*tau;
/*
        printf("a[4] = %f\n",a[4]);
        printf("a[3] = %f\n",a[3]);
        printf("a[2] = %f\n",a[2]);
        printf("a[1] = %f\n",a[1]);
        printf("a[0] = %f\n",a[0]);*/

        x[0] = Aux.get(i,1)-u;
        x[1] = Aux.get(i,2)-u;
        x[2] = Aux.get(i,3)-u;
        x[3] = Aux.get(i,4)-u;
        

    ///////////////// find the first root by Laguerre's method ///////////////////////

        //double x1 = Laguerres(a,4,x[0]);
        double x1 = BirgeVieta(a,4,x[0]);
        double miu1 = x1+u;
        //printf("miu1 = %f\n",miu1);

    ///////////////// find the second root by Laguerre's method //////////////////////
        
        b[3] = a[4];
        b[2] = a[3] + x1*b[3];
        b[1] = a[2] + x1*b[2];
        b[0] = a[1] + x1*b[1];
        //double x2 = Laguerres(b,3,x[1]);
        double x2 = BirgeVieta(b,3,x[1]);
        double miu2 = x2 + u;
        //printf("miu2 = %f\n",miu2);

    ///////////////// find the last two roots by quadratic formula ///////////////////////

        c[2] = b[3];
        c[1] = b[2] + x2*c[2];
        c[0] = b[1] + x2*c[1];
        double x3 = (-c[1]+sqrt(pow(c[1],2)-4.0*c[0]*c[2]))/2.0/c[2];
        double x4 = (-c[1]-sqrt(pow(c[1],2)-4.0*c[0]*c[2]))/2.0/c[2];
        double miu3 = x3 + u;
        double miu4 = x4 + u;
        //printf("miu3 = %f\n",miu3);
        //printf("miu4 = %f\n",miu4);

   ///////////////////////////// find the weights //////////////////////////////////

        y[0]= (m3-m2*miu2-m2*miu3-m2*miu4+m1*miu2*miu3+m1*miu2*miu4+m1*miu3*miu4-rho*miu2*miu3*miu4)
              /((miu1-miu2)*(miu1-miu3)*(miu1-miu4));
        y[1]= -(m3-m2*miu1-m2*miu3-m2*miu4+m1*miu1*miu3+m1*miu1*miu4+m1*miu3*miu4-rho*miu1*miu3*miu4)
              /((miu1-miu2)*(miu2-miu3)*(miu2-miu4));
        y[2]= (m3-m2*miu1-m2*miu2-m2*miu4+m1*miu1*miu2+m1*miu1*miu4+m1*miu2*miu4-rho*miu1*miu2*miu4)
              /((miu1-miu3)*(miu2-miu3)*(miu3-miu4));
        y[3]= -(m3-m2*miu1-m2*miu2-m2*miu3+m1*miu1*miu2+m1*miu1*miu3+m1*miu2*miu3-rho*miu1*miu2*miu3)
              /((miu1-miu4)*(miu2-miu4)*(miu3-miu4));

        
        auxvals.set(i,1,miu1);
        auxvals.set(i,2,miu2);
        auxvals.set(i,3,miu3);
        auxvals.set(i,4,miu4);
        auxvals.set(i,5,y[0]);
        auxvals.set(i,6,y[1]);
        auxvals.set(i,7,y[2]);
        auxvals.set(i,8,y[3]);
        
        }
      
     break;
    }
     
  }

}

//solve polynomials in the form a[n]*x^{n}+a[n-1]*x^{n-1}+...+a[1]*x+a[0] 

double Laguerres(double a[], int n, double z)
  {              
       const int MaxIters = 1000;
       const double TOL = 1.0e-12;
       int mstop = 0;
       int NumIters = 0;
       double d = 0.0;
       double G = 0;
       double H = 0;
       double disc = 0.5;
       

       while (mstop == 0)
         {  
            
            double p = a[n];
            double q = 0.0;
            double r = 0.0;

                   
            for (int i=n-1;i>=0;i--)
              {
                 r = r*z+q;
                 q = q*z+p;
                 p = p*z+a[i];
                 
              }

            G = q/p;
            H = pow(q/p,2)-2.0*r/p;
            disc=(n-1)*(n*H-pow(G,2));


            if (disc>1.0e-12)
              {
                  if ( G > 0 )
                    {
                       d = n/(G+sqrt(disc));
                    }
                  else
                    {
                       d = n/(G-sqrt(disc));    
                    } 
                
                  z = z - d;
                  NumIters = NumIters +1;

                  if (abs(d)<=TOL || NumIters == MaxIters)
                     {
                          mstop = 1;
                     }
              }
            else
               {
                  printf("ERROR in Laguerres: Not all roots are real\n");
                  printf(" disc is %f\n", disc);
                  exit(1);
               }

         }

   return z;
  
}

double BirgeVieta(double a[], int n, double z)
{
       const int MaxIters = 1000;
       const double TOL = 1.0e-12;
       int mstop = 0;
       int NumIters = 0;
       double* r;
       r = new double[n];

     while ( mstop == 0)
       {
         double b = 1.0;
         double c = 1.0;

         for (int i=0;i<n;i++)
          {
              r[i] = a[n-1-i]/a[n];
          }
        for (int i=0;i<n-1;i++)
          {
              b = r[i] + z*b;
              c = b + z*c;
          }
       b = r[n-1]+z*b;

       if ( abs(c) >= TOL )
          {
                z=z-b/c;
          }
       else
          {
            printf(" not all roots are real\n");
            exit(1);
          }
 

       if ( abs(b)<=TOL || NumIters == MaxIters)
          { 
              mstop = 1;
          }
       }
   return z;
}

          

   

