#include "tensors.h"
#include "assert.h"
#include <vector>

using namespace std;

extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

// Update the solution using the constructed Lstar
void UpdateSoln_Interface_ExplicitPart(int run,double alpha1,double alpha2,double beta,double dt,
        const dTensor2& node,const dTensorBC3& aux,
        const dTensorBC3& qstar,const dTensor4& qIstar, const dTensor4 LstarI, const dTensorBC3& Lstar,
        const dTensor4& Implicit,
        dTensorBC3& qnew, dTensor4& qInew,const dTensor4& auxI,
        const dTensorBC1& global2interf,
        const dTensor1& interf2global,
        const dTensor2& dxi)
{

  /*
  int dim=4;
  vector<double> a(dim*dim);
  vector<double> b(dim);
  vector<int> ipiv(dim);
  clock_t t1 = clock();
  srand(1);              // seed the random # generator with a known value
  double maxr=(double)RAND_MAX;
  for(int r=0; r < dim; r++) {  // set a to a random matrix, i to the identity
    for(int c=0; c < dim; c++) {
      a[r + c*dim] = rand()/maxr;
    }
    b[r] = rand()/maxr;
  }
  vector<double> b1(b);
  int info;
  cout << "matrices allocated and initialised " << double(clock() - t1)<< endl;
  clock_t c2 = clock();
	int one = 1;
  dgesv_(&dim, &one, &*a.begin(), &dim, &*ipiv.begin(), &*b.begin(), &dim, &info);
  clock_t c3 = clock();
  cout << "dgesv is over for " << double(c3 - c2) << endl;
  cout << "info is " << info << endl;
  */
    const int     mx = qnew.getsize(1);
    const int   meqn = qnew.getsize(2);
    const int   kmax = qnew.getsize(3);
    const int   maux = aux.getsize(2);
    const int mbc = qnew.getmbc();
    int dim= 6*kmax;
      vector<double> a(dim*dim);
      vector<double> a1(dim*dim);
      vector<double> b(dim);
      vector<int> ipiv(dim);
      vector<double> rhs(dim);


           for(int r=0; r < dim; r++) {  // set a to a random matrix, i to the identity
             for(int c=0; c < dim; c++) {
                a[r+c*dim]=0.0;
                a1[r+c*dim]=0.0;
             }
             b[r]=0.0;
             rhs[r]=0.0;
             ipiv[r]=0;
           }  

//printf("THIS SPOT2! %e \n",global2interf.get(458));
    
#pragma omp parallel for private(a,a1,b,ipiv)
    for (int j=1; j<=(mx+mbc-2); j++)
    {
    //printf("HERE!\n");
    //printf("HERE! %d\n",j);
    //printf("THIS SPOT2! %d \n",j-1);
     
        //printf("THIS SPOT0! %d \n",j-1);
        if( abs(global2interf.get(j-1))<=1.0e-1 && abs(global2interf.get(j))<=1.0e-1 && abs(global2interf.get(j+1))<=1.0e-1)
        {
         for (int m=1; m<=meqn; m++)        
         for (int k=1; k<=kmax; k++)
         {     
          double tmp = alpha1*qstar.get(j,m,k) + alpha2*qnew.get(j,m,k)
            + beta*dt*Lstar.get(j,m,k);
          //printf("HERE! %e %e %e \n",qstar.get(j,m,k),qnew.get(j,m,k),Lstar.get(j,m,k));
          qnew.set(j,m,k, tmp );   
         }
          /*double tmp=0.0;
          if(run==1)
          {    tmp=qstar.get(j,m,k)+0.5*dt*Lstar.get(j,m,k);}
          else
          {    tmp=qnew.get(j,m,k)+dt*Lstar.get(j,m,k);}*/          
        }
 

//printf("THIS SPOT1! %d %d \n",j,(mx+mbc-2));    }
      }
//printf("THIS SPOT3! \n");
    // Optional call to modify updated solution
    void AfterUpdateSoln(const dTensor2& node,
            const dTensorBC3& aux,
            dTensorBC3& q,
            double dt,
            double beta);
    AfterUpdateSoln(node,aux,qnew,dt,beta); 
}

/*
// Update the solution using the constructed Lstar
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensor2& node,const dTensorBC3& aux,
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& q1, dTensorBC3& q2)
{
    const int     mx = q1.getsize(1);
    const int   meqn = q1.getsize(2);
    const int   kmax = q1.getsize(3);
    const int   maux = aux.getsize(2);
    const int mbc = q1.getmbc();

// TODO - do we want to replace these with vget? -DS 06/07/2014
#pragma omp parallel for
    for (int j=(2-mbc); j<=(mx+mbc-1); j++)
    for (int m=1; m<=meqn; m++)        
    for (int k=1; k<=kmax; k++)
    {

        double s1 = q1.get(j,m,k);
        double s3 = qold.get(j,m,k);

        // update q2
        double s2 = q2.get(j,m,k) + delta*s1;
        q2.set(j,m,k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.get(j,m,k);
        q1.set(j,m,k, tmp );
    }

    // Optional call to modify updated solution
    void AfterUpdateSoln(const dTensor2& node,
            const dTensorBC3& aux,
            dTensorBC3& q,
            double dt,
            double beta);
    AfterUpdateSoln(node,aux,q1,dt,beta); 
}*/
