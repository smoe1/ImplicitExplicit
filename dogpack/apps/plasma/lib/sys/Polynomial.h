#ifndef Polynomial_h
#define Polynomial_h
#include "assert.h"
#include "debug.h"

// This class might cause executable bloat
// but inline might be needed for performance.
//
// To enforce positivity in 10-moment model we
// really only need polynomials up to 6th order.
//
#define MAX_ORDER 12
class Polynomial
{
 private:
  double coef[MAX_ORDER];
  // to be really correct the order of the zero
  // polynomial should be -infinity, but 
  // I pretend that it is zero.
  int order;
 private:
  double lead_coef() const
  {
    return coef[order];
  }
  void setzero()
  {
    order=0;
    coef[0]=0.;
  }
  bool vanishing_lead_coef() const
  {
    const double coef_epsilon=1e-12;
    return (lead_coef() < coef_epsilon)
       && (lead_coef() > -coef_epsilon);
  }
  void fix_lead_coef()
  {
    while(order && vanishing_lead_coef())
    {
      order--;
    }
    if(order==0 && vanishing_lead_coef())
    {
      setzero();
    }
  }
 public:
  Polynomial():order(0){coef[0]=1.;}
  Polynomial(const Polynomial& in)
  {
    order = in.order;
    for(int i=0;i<=order;i++) coef[i]=in.coef[i];
  }
  Polynomial& operator=(const Polynomial& in)
  {
    order = in.order;
    for(int i=0;i<=order;i++) coef[i]=in.coef[i];
    return *this;
  }
  Polynomial(double* coef_in, int order_in)
  {
    //assert_le(order,MAX_ORDER);
    //assert_le(0,order);
    order=order_in;
    for(int i=0;i<=order;i++)
    {
      coef[i]=coef_in[i];
    }
    //int i;
    //for(i=0;coef[i]=coef_in[i];i++);
    //order=i-1;
    //if(order<0.) order=0.;
  }
  Polynomial(double c0, double c1)
  {
    //dprint(c0);
    //dprint(c1);
    //assert_ne(c1,0.);
    order = 1;
    coef[0]=c0;
    coef[1]=c1;
  }
  const double get_coef(int i) const
  {
    return coef[i];
  }
  const double* get_coef() const
  {
    return coef;
  }
  void make_monic()
  {
    assert(!iszero());
    double d_inv = 1./lead_coef();
    for(int i=0;i<=order;i++) coef[i]*=d_inv;
  }
  bool iszero() const
  {
    assert_ge(order,0);
    // printf("order==%d\n",order);
    if(order==0 && vanishing_lead_coef()) return true;
    //if(order==0 && coef[0]==0.) return true;
    return false;
  }
  int get_degree() const
  {
    if(order!=0) return order;
    // poor man's substitute for -infinity
    if(!coef[0]) return -1;
    return 0;
  }
  void set_roots(double* roots, int num_roots)
  {
    // choose a monic polynomial.
    order=0;
    coef[0]=1.;
    for(int i=0;i<num_roots;i++)
    {
      multiply_by_linear(roots[i]);
    }
  }
  // multiply by (x-root)
  void multiply_by_linear(double root)
  {
    coef[order+1]=coef[order];
    for(int i=order;i>0;i--)
    {
      coef[i]=coef[i-1]-root*coef[i];
    }
    coef[0]*=-root;
    order++;
  }
  // could define e.g. operator+= instead of these
  void add(double c)
  {
    coef[0]+=c;
  }
  void add(Polynomial const& x)
  {
    if(order>=x.order)
    {
      for(int i=0;i<=x.order;i++)
      {
        coef[i]+=x.coef[i];
      }
    }
    else
    {
      int i;
      for(i=0;i<=order;i++)
      {
        coef[i]+=x.coef[i];
      }
      for(;i<=x.order;i++)
      {
        coef[i]=x.coef[i];
      }
      order = x.order;
    }
  }
  void subtract(double c)
  {
    coef[0]-=c;
  }
  void subtract(Polynomial const& x)
  {
    if(order>=x.order)
    {
      for(int i=0;i<=x.order;i++)
      {
        coef[i]-=x.coef[i];
      }
    }
    else
    {
      int i;
      for(i=0;i<=order;i++)
      {
        coef[i]-=x.coef[i];
      }
      for(;i<=x.order;i++)
      {
        coef[i]=-x.coef[i];
      }
      order = x.order;
    }
  }
  // multiplication
  Polynomial& operator*=(double c)
  {
    if(c==0.)
    {
      coef[0]=0;
      order=0;
    }
    else
    {
      for(int i=0;i<=order;i++)
      {
        coef[i]*=c;
      }
    }
    return *this;
  }
  Polynomial operator*(double c) const
  {
    Polynomial out(*this);
    out*=c;
    return out;
  }
  Polynomial operator*(const Polynomial& x) const
  {
    Polynomial z(*this);
    //z *= x;
    if(z.iszero() || x.iszero())
    {
      z.setzero();
      return z;
    }
    z.order = order+x.order;
    assert_le(z.order,MAX_ORDER);
    for(int k=0;k<=z.order;k++) z.coef[k]=0.;

    for(int i=0;i<=x.order;i++)
    for(int j=0;j<=order;j++)
    {
      z.coef[i+j]+=x.coef[i]*coef[j];
    }
    return z;
  }
  Polynomial operator-()
  {
    Polynomial out;
    out.order=order;
    for(int i=0;i<=out.order;i++)
    {
      out.coef[i]=-coef[i];
    }
    return out;
  }
  Polynomial& operator/(double s)
  {
    assert_ne(s,0.);
    for(int i=0;i<=order;i++) coef[i]/=s;
  }
  Polynomial& operator*=(Polynomial const& x)
  {
    *this = (*this) * x;
    return *this;
  }
  double eval(double x) const
  {
    // evaluate using Horner's method
    double out=coef[order];
    for(int i=order-1;i>=0;i--)
    {
      out *= x;
      out += coef[i];
    }
    return out;
  }
  // divide by x-c (using synthetic division) and return the remainder
  // (which is simply the value of the polynomial at c)
  double get_remainder_set_quotient(double c)
  {
    // In this case we simply end up shifting the coefficients
    // and returning the constant term
    //assert_ne(c,0.);
    if(order==0)
    {
      double out=coef[0];
      coef[0]=0.;
      return out;
    }
    double out=coef[order--];
    for(int i=order;i>=0;i--)
    {
      double prod = out*c;
      double sum = prod + coef[i];
      coef[i] = out;
      out = sum;
    }
    return out;
  }
  // divide: discard remainder and set to quotient.
  void set_quotient(const Polynomial& divisor)
  {
    // This implementation is slightly wasteful.
    Polynomial remainder = *this;
    *this = remainder.get_quotient_set_remainder(divisor);
  }
  Polynomial get_remainder_set_quotient(const Polynomial& divisor)
  {
    Polynomial remainder = *this;
    *this = remainder.get_quotient_set_remainder(divisor);
    return remainder;
  }
  // divide: return quotient and replace with remainder
  Polynomial get_quotient_set_remainder(const Polynomial& divisor)
  {
    //assert(!divisor.iszero());
    double d = divisor.lead_coef();
    assert_ne(d,0.);
    if(divisor.order==0)
    {
      Polynomial quotient=*this/d;
      setzero();
      return quotient;
    }
    int k = order - divisor.order;
    Polynomial quotient;
    quotient.order=k;
    for(;k>=0;k--)
    {
      // the coefficient of the quotient
      double q_k = lead_coef()/d;
      quotient.coef[k] = q_k;
      // remainder := remainder - q_k*x^k*divisor
      order--;
      for(int m=0;m<divisor.order;m++)
      {
        coef[k+m] -= q_k*divisor.coef[m];
      }
    }
    // decrement order of out until its
    // leading coefficient is nonzero
    // or it is the zero polynomial.
    fix_lead_coef();
    return quotient;
  }
  // divide, returning the remainder and discarding the quotient.
  Polynomial get_remainder(const Polynomial& divisor) const
  {
    //assert(!divisor.iszero());
    double d = divisor.lead_coef();
    assert_ne(d,0.);
    Polynomial out(*this);
    if(divisor.order==0)
    {
      out.setzero();
      return out;
    }
    for(int k = order-divisor.order;k>=0;k--)
    {
      // the (unretained) coefficient of the quotient
      double q_k = out.lead_coef()/d;
      // out := out - q_k*x^k*divisor
      out.order--;
      for(int m=0;m<divisor.order;m++)
      {
        out.coef[k+m] -= q_k*divisor.coef[m];
      }
    }
    // decrement order of out until its
    // leading coefficient is nonzero
    // or it is the zero polynomial.
    out.fix_lead_coef();
    return out;
  }
  Polynomial get_derivative() const
  {
    Polynomial out(*this);
    out.differentiate();
    return out;
  }
  void differentiate()
  {
    if(order==0)
    {
      coef[0]=0.;
      return;
    }
    int i=0;
    int j=1;
    for(;j<=order;i++,j++)
    {
      coef[i]=coef[j]*j;
    }
    coef[order]=0.;
    order--;
  }

  //double ridders_root(double low, double hgh);
  void isolate_root(double& low, double& hgh);
  void isolate_and_divide_out_root_minus_x(double& low, double& hgh)
  {
    isolate_root(low, hgh);
    double remainder = get_remainder_set_quotient(low);
    //dprint3(remainder);
    (*this)*=-1;
  }
  double polish_root(double x0) const;
  void polish_roots(double* roots, int num_roots) const;
  double get_smallest_positive_root() const;

  void print_matlab()
  {
    printf("[%24.16f",coef[order]);
    for(int i=order-1;i>=0;i--)
    {
      printf(" %24.16f",coef[i],i);
    }
    printf("]");
  }
  void printx()
  {
    printf("%.1f",coef[0]);
    for(int i=1;i<=order;i++)
    {
      printf("%+.1fx^%d",coef[i],i);
    }
  }
  void print()
  {
    print_matlab();
  }
};

inline Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs)
{
  Polynomial temp = lhs;
  temp.add(rhs);
  return temp;
}

inline Polynomial operator-(const Polynomial& lhs, double rhs)
{
  Polynomial temp = lhs;
  temp.subtract(rhs);
  return temp;
}

inline Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs)
{
  Polynomial temp = lhs;
  temp.subtract(rhs);
  return temp;
}

inline Polynomial operator*(double lhs, const Polynomial& rhs)
{
  return rhs*lhs;
}

#endif
