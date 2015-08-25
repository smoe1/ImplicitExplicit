#include "math.h"
#include "../Polynomial.h"
#include "../JTrootFinder.h"

void test2()
{
  double pc[12] = {4, 3, 2, 0};
  double bc[12] = {-5,1,0};
  Polynomial p(pc,2);
  Polynomial pcopy(p);
  Polynomial b(bc,1);
  p.print();
  printf(" at x=5 equals %.0f\n", p.eval(5));
  double remainder = p.get_remainder_set_quotient(5);
  printf("That is, divided by (x-5) gives remainder %.2f\n", remainder);
  printf("That is, ");
  pcopy.print();
  printf(" equals ");
  b.print();
  printf(" times ");
  // printf(" something ");
  Polynomial q = pcopy.get_quotient_set_remainder(b);
  q.print();
  printf(" plus ");
  //Polynomial r = pcopy.remainder(b);
  //r.print();
  pcopy.print();
  printf("\n");
  printf("Indeed, we get ");
  (b*q+pcopy).print();
  printf("\n");
}

void test1()
{
    double c1[6] = {2,-1,1,0};
    double c2[6] = {3,2,1,0};
    Polynomial p1(c1,2);
    Polynomial p2(c2,2);
    //
    p1.print();
    printf(" + ");
    p2.print();
    printf(" = ");
    (p1+p2).print();
    printf("\n");
    //
    p1.print();
    printf(" * ");
    p2.print();
    printf(" = ");
    Polynomial p3 = p1*p2;
    p3.print();
    printf("\n");

    Polynomial q = p3.get_quotient_set_remainder(p1);
    printf(" p3 should now be zero: ");
    p3.print();
    printf("\nAnd q*p1 should give the original product: ");
    Polynomial p3test = q*p1;
    p3test.print();
    printf("\n");
    //
    printf("derivative of ");
    p2.print();
    printf(" is ");
    p2.differentiate();
    p2.print();
    printf("\n");
}

bool sign_change(double old_val, double new_val)
{
  if(old_val <= 0. && new_val > 0.) return true;
  if(old_val >= 0. && new_val < 0.) return true;
  return false;
}

int count_sign_changes(double x, Polynomial* sturm, int num_sturm)
{
  int count=0;
  double val = sturm[0].eval(x);
  for(int k=1;k<num_sturm;k++)
  {
    double val_old = val;
    val = sturm[k].eval(x);
    if(sign_change(val_old,val)) count++;
  }
  return count;
}

int count_roots_in_interval(double a, double b,
  Polynomial* sturm, int num_sturm)
{
  int ca = count_sign_changes(a, sturm, num_sturm);
  int cb = count_sign_changes(b, sturm, num_sturm);
  dprint(ca);
  dprint(cb);
  return ca-cb;
}

void print_num_roots(double a,double b,Polynomial* sturm,int num_sturm)
{
  printf("The number of roots in [%f,%f) is %d\n", a,b,
    count_roots_in_interval(a, b, sturm, num_sturm));
}

void test_sturm()
{
  double r[6] = {3, 3, 4, 4, 5, 0};
  int num_roots = 5;
  // build sturm sequence
  //
  Polynomial sturm[8];
  sturm[0].set_roots(r,num_roots);
  sturm[1] = sturm[0].get_derivative();
  // number of polynomials in sturm sequence.
  printf("n=%d\n",0);
  sturm[0].print(); printf("\n");
  int n=1;
  if(!sturm[1].iszero()) while(1)
  {
    printf("n=%d\n",n);
    sturm[n].print(); printf("\n");
    dprint(sturm[n].iszero());
    sturm[n+1] = -sturm[n-1].get_remainder(sturm[n]);
    n++;
    if(sturm[n].iszero()) break;
    // does not seem to help with stability
    sturm[n].make_monic();
  }
  int num_sturm = n;
  // we now verify that sturm[num_sturm-1] is
  // a constant polynomial; else it is the gcd of
  // the sturm sequence and we should divide it out
  // of all the polynomials in the sequence
  // (which does not get rid of any roots).
  // Note that dividing out the gcd does not affect
  // the number of sign changes except when evaluating
  // precisely at a root of the gcd.
  Polynomial gcd = sturm[num_sturm-1];
  printf("gcd(p,p') = ");
  gcd.print(); printf("\n");
  int gcd_degree = gcd.get_degree();
  #if 1
  if(gcd_degree>0)
  {
    for(int n=0;n<num_sturm;n++)
    {
      sturm[n].print(); printf("\n");
      printf("remainder = ");
      Polynomial remainder = sturm[n].get_remainder_set_quotient(gcd);
      remainder.print(); printf("\n");
      //assert(remainder.iszero());
      sturm[n].print(); printf("\n");
    }
  }
  assert_eq(sturm[num_sturm-1].get_degree(),0);
  #endif
  // (if we cared about repeated zeros we could
  // compute the sturm sequence for the gcd (recursively).)

  // Use sturm sequence to estimate number of zeros in interval
  //
  #if 0
  double a=3.;
  double b=7.0000;
  int ca=0;
  int cb=0;
  double val_a=sturm[0].eval(a);
  double val_b=sturm[0].eval(b);
  for(int k=1;k<num_sturm;k++)
  {
    double val_a_old = val_a;
    double val_b_old = val_b;
    val_a = sturm[k].eval(a);
    val_b = sturm[k].eval(b);
    if(sign_change(val_a_old,val_a)) ca++;
    if(sign_change(val_b_old,val_b)) cb++;
  }
  #endif
  printf("roots: ");
  for(int i=0;i<num_roots;i++) printf("%f, ",r[i]);
  printf("\n");
  //printf("The number of roots in [%f,%f) is %d-%d = %d\n",
  //  a,b,ca,cb,ca-cb);
  print_num_roots(2.9,5.1,sturm,num_sturm);
}

void test_repeated_real_roots()
{
  //const int ORDER = 8;
  //double r[ORDER] = {1, 2, 2, 2, 3, 3, 7, 7};
  const int ORDER = 6;
  double r[ORDER] = {1, 2, 2, 2, 3, 3};
  double zeror[ORDER], zeroi[ORDER];
  int num_roots = ORDER;
  int info[16];
  double c[ORDER];
  //
  Polynomial poly;
  poly.set_roots(r,num_roots);
  const double* coef = poly.get_coef();
  const int degree = poly.get_degree();
  assert_eq(degree,num_roots);
  for(int i=0;i<=degree;i++) c[i] = coef[degree-i];
  JTrootFinder jtRootFinder;
  int num_roots2 = jtRootFinder.rpoly(c, num_roots, zeror, zeroi, info);
  assert_eq(num_roots2,num_roots);
  //
  printf("roots: ");
  for(int i=0;i<num_roots;i++) printf("%24.16e, ",zeror[i]);
  printf("\n");
  printf("values at roots: ");
  for(int i=0;i<num_roots;i++) printf("%24.16e, ",poly.eval(zeror[i]));
  printf("\n");
  //
  #if 0
  printf("roots: ");
  for(int i=0;i<num_roots;i++) zeror[i] = poly.polish_root(zeror[i]);
  for(int i=0;i<num_roots;i++) printf("%24.16e, ",zeror[i]);
  printf("\n");
  printf("values at roots: ");
  for(int i=0;i<num_roots;i++) printf("%24.16e, ",poly.eval(zeror[i]));
  printf("\n");
  #endif
}

void test_roots()
{
  double c[5]={-3,1};
  Polynomial p(c,1);
  p.multiply_by_linear(1.);
  p.print();
  printf("\n");
}

void test_isolate_root()
{
  const int ORDER = 5;
  double r[ORDER] = {1.5, 1.7, 2.29, 2.3, 3.3};
  //double r[ORDER] = {1.5, 1.7, 2.29, 2.3, 2.31};
  //double r[ORDER] = {1.5, 1.7, 2.299999, 2.3, 2.300001};
  //double r[ORDER] = {1.5, 1.7, 2.2999999, 2.3, 2.3000001};
  double zeror[ORDER], zeroi[ORDER];
  int num_roots = ORDER;
  int info[16];
  double c[ORDER];
  //
  Polynomial poly;
  poly.set_roots(r,num_roots);
  double root = poly.isolate_root(0.,4.7);
  //
  printf("%24.16e, ",root);
  printf("\n");
  printf("value at root: ");
  printf("%24.16e, ",poly.eval(root));
  printf("\n");
}

int 
/*
gsl_poly_solve_cubic (double a, double b, double c, 
                      double *x0, double *x1, double *x2);
*/
gsl_poly_solve_cubic (const double* p, double *roots);
void test_solve_cubic()
{
  const int ORDER = 3;
  // double r[ORDER] = {1.5, 1.5, 1.7}; // gives 1 root
  double r[ORDER] = {-1.3, -0.299999, -0.299999};
  double zeros[ORDER];
  //
  Polynomial poly;
  poly.set_roots(r,ORDER);
  const double* coef = poly.get_coef();
  double rt[ORDER];
  /*
  int numroots = gsl_poly_solve_cubic(coef[2],coef[1],coef[0],
    &rt[0],&rt[1],&rt[2]);
  */
  int numroots = gsl_poly_solve_cubic(coef, rt);
  //
  dprint(numroots);
  printf("roots: \n");
  dprint(rt[0]);
  dprint(rt[1]);
  dprint(rt[2]);
  printf("values at root: \n");
  dprint(poly.eval(rt[0]));
  dprint(poly.eval(rt[1]));
  dprint(poly.eval(rt[2]));
}

int main()
{
  //test1();
  //test2();
  //test_sturm();
  //test_repeated_real_roots();
  //test_isolate_root();
  test_solve_cubic();
  //test_roots();
}
