#include "Source2Fluid.h"
#include <ctime>
#include <cstdlib>

#include <stdlib.h> // for rand, srand
#include <time.h> // to seed

static void seed_rand()
{
  srand((unsigned)time(NULL)); 
}

static double rand_double()
{
  double out;
  // produces a number in [0, 1)
  double scale_inv = 1./(RAND_MAX+1.);
  //out = double(rand())/(RAND_MAX+1.);
  out = double(rand())*scale_inv;
  return out;
}

void test_set_B_orth_system()
{
  seed_rand();
  const double B1 = rand_double();
  const double B2 = rand_double();
  const double B3 = rand_double();
  const double Bmag_inv = 1./sqrt(B1*B1+B2*B2+B3*B3);
  a[4][4];
  set_B_orth_system(double a[4][4],B1,B2,B3,Bmag_inv);
}

int main()
{
  test_set_B_orth_system();
}
