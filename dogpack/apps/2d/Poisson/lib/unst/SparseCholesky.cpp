// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (SparseCholesky.cpp)
// --------------------------------------------------------------------------

#include "SparseCholesky.h"
#include <cstdlib>
#include <cstdio>
#include "assert.h"
#include "tensors.h"

SparseCholesky::SparseCholesky(int size_in)
{
  size = size_in;
  num_nonzeros = 0;

  Rnz_row = new int[size+1];
  Rind    = new int*[size+1];
  Rval    = new double*[size+1];

  Lnz_row = new int[size+1];
  Lind    = new int*[size+1];
  Lval    = new double*[size+1];

  P = new int[size+1];
}

SparseCholesky::SparseCholesky(const SparseCholesky& amat)
{
  size = amat.size;
  num_nonzeros = amat.num_nonzeros;

  for (int i=1; i<=size; i++)
    {
      Rnz_row[i] = amat.Rnz_row[i];
      Lnz_row[i] = amat.Lnz_row[i];
      P[i] = amat.P[i];
    }

  for (int i=1; i<=size; i++)
    for (int j=1; j<=Rnz_row[i]; j++)
      {
	Rind[i][j] = amat.Rind[i][j];
	Rval[i][j] = amat.Rval[i][j];
      }

  for (int i=1; i<=size; i++)
    for (int j=1; j<=Lnz_row[i]; j++)
      {
	Lind[i][j] = amat.Lind[i][j];
	Lval[i][j] = amat.Lval[i][j];
      }
}

SparseCholesky::~SparseCholesky()
{
  delete[] Rnz_row;
  Rnz_row = NULL;

  delete[] Rind;
  Rind = NULL;

  delete[] Rval;
  Rval = NULL;

  delete[] Lnz_row;
  Lnz_row = NULL;

  delete[] Lind;
  Lind = NULL;

  delete[] Lval;
  Lval = NULL;

  delete[] P;
  P = NULL;
}


const double& SparseCholesky::get_Rval(int i, int j) const
// Returns jth non-zero entry on ith row
{
  return Rval[i][j];
}

const int& SparseCholesky::get_Rind(int i, int j) const
// Returns column index of jth non-zero entry on ith row
{
  return Rind[i][j];
}

const double& SparseCholesky::get_Lval(int i, int j) const
// Returns jth non-zero entry on ith row
{
  return Lval[i][j];
}

const int& SparseCholesky::get_Lind(int i, int j) const
// Returns column index of jth non-zero entry on ith row
{
  return Lind[i][j];
}

const int& SparseCholesky::get_Pval(int i) const
// Returns ith entry of permutation array
{
  return P[i];
}

void SparseCholesky::read(const char* inputdir)
// Read in sparse matrix
//
{
  char fname[1024];

  // input: number of non-zeros in R
  snprintf(fname,1024,"%s/matlab_output/Rnz.dat",inputdir);
  FILE* read1 = fopen(fname,"r");
  if (read1==NULL)
    { filenotfound(fname); }
  num_nonzeros = 0;
  Rnz_row[0]=0;

  for (int i=1; i<=size; i++)
    {
      int tmp_int;

      fscanf(read1,"%i",&tmp_int);
      Rnz_row[i] = tmp_int;  

      num_nonzeros = num_nonzeros + tmp_int;

      Rind[i] = new int[tmp_int+1];
      Rval[i] = new double[tmp_int+1];
    }
  fclose(read1);

  // input: get all non-zeros in R (index and value)
  int* counter = new int[size+1];
  for (int i=0; i<=size; i++)
    {  counter[i] = 0;  }
  snprintf(fname,1024,"%s/matlab_output/R.dat",inputdir);
  FILE* read2 = fopen(fname,"r");
  if (read2==NULL)
    { filenotfound(fname); }
  for (int m=1; m<=num_nonzeros; m++)
    {
      int i,j;
      double tmp_double;

      fscanf(read2,"%i %i %lf",&i,&j,&tmp_double);
      counter[i] = counter[i] + 1;

      Rind[i][counter[i]] = j;
      Rval[i][counter[i]] = tmp_double;
    }
  fclose(read2);

  // input: number of non-zeros in L
  snprintf(fname,1024,"%s/matlab_output/Lnz.dat",inputdir);
  FILE* read3 = fopen(fname,"r");
  if (read3==NULL)
    { filenotfound(fname); }
  int num_nonzeros_tmp = 0;
  Lnz_row[0]=0;
  for (int i=1; i<=size; i++)
    {
      int tmp_int;

      fscanf(read2,"%i",&tmp_int);
      Lnz_row[i] = tmp_int;  

      num_nonzeros_tmp = num_nonzeros_tmp + tmp_int;

      Lind[i] = new int[tmp_int+1];
      Lval[i] = new double[tmp_int+1];
    }
  fclose(read3);
  assert_eq(num_nonzeros,num_nonzeros_tmp);

  // input: get all non-zeros in L (index and value)
  for (int i=0; i<=size; i++)
    {  counter[i] = 0;  }
  snprintf(fname,1024,"%s/matlab_output/L.dat",inputdir);
  FILE* read4 = fopen(fname,"r");
  if (read4==NULL)
    { filenotfound(fname); }
  for (int m=1; m<=num_nonzeros; m++)
    {
      int i,j;
      double tmp_double;
      fscanf(read4,"%i %i %lf",&i,&j,&tmp_double);

      counter[i] = counter[i] + 1;

      Lind[i][counter[i]] = j;
      Lval[i][counter[i]] = tmp_double;
    }
  fclose(read4);

  // input: get all permutation values
  snprintf(fname,1024,"%s/matlab_output/P.dat",inputdir);
  FILE* read5 = fopen(fname,"r");
  if (read5==NULL)
    { filenotfound(fname); }
  P[0]=0;
  for (int m=1; m<=size; m++)
    {
      int tmp_int;
      fscanf(read4,"%i",&tmp_int);

      P[m] = tmp_int;
    }
  fclose(read5);
}

void SparseCholesky::write(const char* outputdir)
// Write out sparse matrix
{
  char fname[1024];

  // output: number of non-zeros in R
  snprintf(fname,1024,"%s/matlab_output/Rnz_out.dat",outputdir);
  FILE* write1 = fopen(fname,"w");
  for (int i=1; i<=size; i++)
    {
      fprintf(write1,"%8i\n",Rnz_row[i]);
    }
  fclose(write1);

  // output: write all non-zeros in R (index and value)
  snprintf(fname,1024,"%s/matlab_output/R_out.dat",outputdir);
  FILE* write2 = fopen(fname,"w");
  for (int i=1; i<=size; i++)
    for (int k=1; k<=Rnz_row[i]; k++)
      {
	int j = Rind[i][k];
	
	fprintf(write2,"%8i %8i %32.16e\n",i,j,Rval[i][k]);
      }
  fclose(write2);

  // output: number of non-zeros in L
  snprintf(fname,1024,"%s/matlab_output/Lnz_out.dat",outputdir);
  FILE* write3 = fopen(fname,"w");
  for (int i=1; i<=size; i++)
    {
      fprintf(write3,"%8i\n",Lnz_row[i]);
    }
  fclose(write3);

  // output: write all non-zeros in L (index and value)
  snprintf(fname,1024,"%s/matlab_output/L_out.dat",outputdir);
  FILE* write4 = fopen(fname,"w");
  for (int i=1; i<=size; i++)
    for (int k=1; k<=Lnz_row[i]; k++)
      {
	int j = Lind[i][k];
	fprintf(write4,"%8i %8i %32.16e\n",i,j,Lval[i][k]);
      }
  fclose(write4);

  // output: write all permutation values
  snprintf(fname,1024,"%s/matlab_output/P_out.dat",outputdir);
  FILE* write5 = fopen(fname,"w");
  P[0]=0;
  for (int m=1; m<=size; m++)
    {
      fprintf(write5,"%8i\n",P[m]);
    }
  fclose(write5);
}

const int& SparseCholesky::get_size() const
// Returns matrix size (i.e., number of rows = number of columns)
{
  return size;
}

const int& SparseCholesky::get_num_nonzeros() const
// Returns total number of nonzero entries 
//   0 <= num_nonzeros <= size*size
{
  return num_nonzeros;
}

const int& SparseCholesky::get_Rnz_row(int i) const
{
  return Rnz_row[i];
}

const int& SparseCholesky::get_Lnz_row(int i) const
{
  return Lnz_row[i];
}

const void SparseCholesky::ForwardSubs(const dTensor1& b, dTensor1& y) const
// Forward substitution using L = R^T
{
  for (int i=1; i<=size; i++)
    {
      int kend = Lnz_row[i];
      double denom = Lval[i][kend]; 
      double numer = b.get(P[i]);

      //assert_eq(Lind[i][kend],i);

      for (int k=1; k<=(kend-1); k++)
        {
	  int j      = Lind[i][k];
	  double val = Lval[i][k];
	  numer = numer - val*y.get(j);
        }

      y.set(i, numer/denom );
    }
}

const void SparseCholesky::BackwardSubs(const dTensor1& y, dTensor1& x) const
// Backward substitution using R
{
  dTensor1 xtmp(size);

  for (int i=size; i>=1; i--)
    {
      double denom = Rval[i][1]; 
      double numer = y.get(i);

      //assert_eq(Rind[i][1],i);

      for (int k=2; k<=Rnz_row[i]; k++)
        {
	  int j      = Rind[i][k];
	  double val = Rval[i][k];
	  numer = numer - val*xtmp.get(j);
        }

      xtmp.set(i, numer/denom );
    }

  for (int i=1; i<=size; i++)
    {
      x.set(P[i], xtmp.get(i) );
    }

}

void SparseCholesky::filenotfound(char* fname)
{
  printf("\n");
  printf("  Error in SparseCholesky.cpp.\n");
  printf("  File not found: %s.\n",fname);
  printf("\n");
  exit(1);
}
