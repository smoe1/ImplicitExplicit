//#define MAX_DEBUG_LEVEL 3 //comment out to turn off verbose debug in this file
#include "Polynomial.h"
#include "JTrootFinder.h"
#include "math.h"
#include "debug.h"
#include "float.h" // DBL_MAX

#define GSL_DBL_EPSILON 2.2204460492503131e-16

static inline double Min(double a, double b)
{ return a<b?a:b; }

class PointVal
{
  double _x;
  double _v; // eval(x)
  double _a; // fabs(v)
 public:
  PointVal(){}
  PointVal(const PointVal& in):
   _x(in._x),
   _v(in._v),
   _a(in._a)
  {}
  PointVal(double xin, const Polynomial* p)
  {
    eval(xin, p);
  }
  void eval(double xin, const Polynomial* p)
  {
    _x=xin;
    _v=p->eval(_x);
    _a=fabs(_v);
  }
  double x() const {return _x;}
  double v() const {return _v;}
  double a() const {return _a;}
};

class RootBracket
{
  PointVal _low; // maintain low().x() > 0;
  PointVal _hgh; // maintain hgh().x() < 0;
 public:
  const PointVal& low() const {return _low;}
  const PointVal& hgh() const {return _hgh;}
  void low(const PointVal& in) {_low=in;}
  void hgh(const PointVal& in) {_hgh=in;}
  bool not_in_interval(double x)
  {
    if(x > hgh().x() || x < low().x()) return true;
    return false;
  }
};

class RootBracketState : public RootBracket
{
 private:
  // the most recently evaluated point
  PointVal _pv;
  int _num_evals;
  // the most recently updated end
  bool _using_hgh;
  double var_epsilon;
  bool _done;
 public:
  const PointVal& pv() const {return _pv;}
  const PointVal& pv_old() const {
    if(_using_hgh) return low();
    return hgh();
  }
  bool done() const {return _done;}
  bool using_hgh() const {return _using_hgh;}
  void set_to_using_smaller()
  {
    if(hgh().a() < low().a())
    {
      _using_hgh = true;
      _pv=hgh();
    }
    else
    {
      _using_hgh = false;
      _pv=low();
    }
  }
  double x() const {return _using_hgh ? hgh().x() : low().x();}
  RootBracketState(const PointVal& pv_low, const PointVal& pv_hgh):
    _num_evals(2),var_epsilon(1e-16),_done(false)
  {
    low(pv_low);
    hgh(pv_hgh);
    if(low().v()*hgh().v() > 0.)
    {
      dprint(low().x());
      dprint(hgh().x());
      dprint(low().v());
      dprint(hgh().v());
      eprintf("signs agree");
    }
    set_to_using_smaller();
  }
  double get_ratio() const {
    return _using_hgh ? hgh().a()/low().a()
                      : low().a()/hgh().a();
  }
  int num_evals() const {return _num_evals;}
  // evaluate at a point in the interval
  bool eval(double x, Polynomial* p)
  {
    _pv.eval(x,p);
    _num_evals++;
    dprintf3("iter %d: value at x=%24.16e is %24.16e",_num_evals,x,_pv.v());
    if(_num_evals > 30)
    {
      var_epsilon *= 2;
      assert_lt(_num_evals,100); // pure bisection does better than this
    }
    if(_pv.a() < var_epsilon) _done=true;
    //
    // confirm that the point is in the interval
    //assert_gt(x,low().x()); assert_lt(x,hgh().x());
    if(x<=low().x())
    {
      dprintf("x=%24.16e falls below lower bound %24.16e",x,low().x());
      return _done;
    }
    if(x>=hgh().x())
    {
      dprintf("x=%24.16e lands above upper bound %24.16e",x,hgh().x());
      return _done;
    }
    if(_pv.v() > 0.)
    {
      low(_pv);
      _using_hgh = false;
    }
    else
    {
      hgh(_pv);
      _using_hgh = true;
    }
    return _done;
  }
};

class RootBracketClusterState : public RootBracketState
{
 private:
  int _num_roots_in_cluster;
 public:
  int num_roots_in_cluster() const {return _num_roots_in_cluster;}
  void reset_num_roots_in_cluster() {
    _num_roots_in_cluster = 1;
    dprint(_num_roots_in_cluster);
  }
  void increment_num_roots_in_cluster()
  {
    _num_roots_in_cluster +=2;
    dprint(_num_roots_in_cluster);
  }
  void decrement_num_roots_in_cluster()
  {
    assert(_num_roots_in_cluster >= 3);
    _num_roots_in_cluster -=2;
    dprint(_num_roots_in_cluster);
  }
  RootBracketClusterState(const PointVal& pv_low, const PointVal& pv_hgh):
    RootBracketState(pv_low,pv_hgh), _num_roots_in_cluster(1) { }
};

// class NewtonState
// {
//   bool using_hgh;
//   // the current point
//   PointVal x0;
//  public:
// };
// 
// class SecantState
// {
//   // the current points
//   PointVal x0;
//   PointVal x1;
// };

// the trick is to set the trigger small enough to guarantee that
// if the assumed number of roots is incorrect the trigger will
// eventually be triggered.
//
inline double get_trigger_ratio_for_newton(int num_roots_in_cluster)
{
  const double threshhold = 0.5; // 0 < 0.5 < 1
  double trigger_ratio = threshhold*exp(-1);
  if(num_roots_in_cluster!=1)
  {
    switch(num_roots_in_cluster)
    {
     case 1:
     case 2:
     case 4:
      eprintf("disallowed");
      break;
     case 3:
      trigger_ratio = threshhold*exp(-3);
      break;
     case 5:
      trigger_ratio = threshhold*exp(-5);
      break;
     default:
      // exp is usually implemented with lookup table nowadays
      trigger_ratio = threshhold*exp(-num_roots_in_cluster);
      break;
    }
  }
  return trigger_ratio;
}

// can we show that the maximal effective trigger for 
// newton is always greater than for secant method?
//
// in the case of the secant method, for the maximal trigger the
// ratio of consecutive residuals eventually alternates between
// being lower and higher than the trigger.
//
inline double get_trigger_ratio_for_secant(int num_roots_in_cluster)
{
  switch(num_roots_in_cluster)
  {
   // I want to avoid increasing num_roots_in_cluster from 1
   // unless it is really called for
   case 1:
    return .4; // 0 < .4 < pow((sqrt(5.)-1.)/2.,1.5) = .4859 = maximal trigger
    break;
  }
  return get_trigger_ratio_for_newton(num_roots_in_cluster);
}

// ratio = v0/v1.
inline double secant_rule(double x0,double x1,double ratio,
  int num_roots_in_cluster)
{
    switch(num_roots_in_cluster)
    {
     case 1:
      break;
     case 2:
      eprintf("disallowed");
      break;
     case 3:
      ratio = cbrt(ratio);
      break;
     default:
      assert_eq(1,num_roots_in_cluster%2);
      ratio = copysign(pow(fabs(ratio),1./num_roots_in_cluster),ratio);
    }
    double denominator = (1.-ratio);
    assert_gt(denominator,0.);
    return (x0-ratio*x1)/denominator;
}

// #include"nrutil.h" 
// #define MAXIT 60 
// #define UNUSED(-1.11e30) 
// double Polynomial::ridders_root(double x1, double x2)
// {
//   double xacc = 1e-16;
//   int j; 
//   double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew; 
//   fl=eval(x1); 
//   fh=eval(x2); 
//   if((fl>0.0&&fh< 0.0)||(fl< 0.0&&fh>0.0)){ 
//     xl=x1; 
//     xh=x2; 
//     ans=UNUSED;
//     
//     for(j=1;j<=MAXIT;j++){ 
//       xm=0.5*(xl+xh); 
//       fm=eval(xm);
//       s=sqrt(fm*fm-fl*fh); 
//       if(s==0.0) return ans; 
//       xnew=xm+(xm-xl)*((fl>=fh?1.0: -1.0)*fm/s);
//       if(fabs(xnew-ans)<=xacc) return ans; 
//       ans=xnew; 
//       fnew=eval(ans);
//       if(fnew==0.0) return ans; 
//       if(SIGN(fm,fnew)!=fm){
//         xl=xm; 
//         fl=fm; 
//         xh=ans; 
//         fh=fnew; 
//       } else if(SIGN(fl,fnew)!=fl){ 
//         xh=ans; 
//         fh=fnew; 
//       } else if(SIGN(fh,fnew)!=fh){ 
//         xl=ans; 
//         fl=fnew; 
//       } else eprintf("never get here."); 
//       if(fabs(xh-xl)<=xacc) return ans; 
//     } 
//     eprintf("zriddr exceed maximum iterations"); 
//   } 
//   else{ 
//     if(fl==0.0) return x1; 
//     if(fh==0.0) return x2; 
//     eprintf("root must be bracketed in zriddr."); 
//   } 
//   return 0.0;
// }

// double Polynomial::my_ridders_root(double low, double hgh)
// {
//  if(0)
//  {
//   dprint(low);
//   dprint(hgh);
//   // ensure starting conditions
//   //
//   PointVal pv_low(low,this);
//   // double p_low = eval(low);
//   if(pv_low.v() < 0.)
//   {
//     Polynomial p_minus = -(*this);
//     return p_minus.isolate_root(low,hgh);
//   }
// 
//   // initialize algorithm
//   //
//   // initialize root bracket.
//   //
//   PointVal pv_hgh(hgh,this);
//   RootBracketState rb(pv_low,pv_hgh);
//   dprint(rb.low().x());
//   dprint(rb.low().v());
//   dprint(rb.hgh().x());
//   dprint(rb.hgh().v());
//   int num_evals = 0;
//  
//  ridders_iteration:
//   const double x1 = rb.low().x();
//   const double x2 = rb.hgh().x();
//   const double f1 = rb.low().v();
//   const double f2 = rb.hgh().v();
//   const double x3 = (x1+x2)/2.; // midpoint
//   if(rb.eval(x3,this)) goto done;
//   const double f3 = rb.pv().v();
//   const double x4 = x3+(x3-x1)*(copysign(f3,f1-f2)/sqrt(pow(f3,2)-f1*f2));
//   if(rb.eval(x4,this)) goto done;
//   goto ridders_iteration;
// 
//  done:
//    dprint(rb.num_evals());
//    return rb.pv().x();
//  }
// }

double isolate_root_using_bisection(double& low, double& hgh, Polynomial *p)
{
  double x0, p0;

  int bisection_iterations=0;
  double p_low = p->eval(low);
  if(p_low < 0.){
    Polynomial minus_p = -*p;
    return isolate_root_using_bisection(low, hgh, &minus_p);
  }
  double p_hgh = p->eval(hgh);
  assert_ge(p_low,0.);
  assert_le(p_hgh,0.);
  // for smaller intervals average is not necessarily
  // strictly between endpoints of interval
  double x_epsilon = 4*GSL_DBL_EPSILON;
  double v_epsilon = 4*GSL_DBL_EPSILON;
 bisect:
  {
    bisection_iterations++;
    assert_lt(bisection_iterations,100);
    x0 = 0.5*(low+hgh);
    //dprint(x0);
    p0 = p->eval(x0);
    //dprint(p0);
    //const double ap0 = fabs(p0);
    //if(ap0<v_epsilon) goto done;
    if(p0<0.){
      hgh = x0;
      //dprintf("low=%24.16e, hgh=%24.16e\n",low,hgh);
      double interval_length = x0-low;
      //dprint(interval_length);
      if(interval_length<x_epsilon) goto done;
    } else {
      low = x0;
      //dprintf("low=%24.16e, hgh=%24.16e\n",low,hgh);
      double interval_length = hgh-x0;
      //dprint(interval_length);
      if(interval_length<x_epsilon) goto done;
      if(p0<v_epsilon) goto done;
    }
    goto bisect;
  }
  done:
    dprint(bisection_iterations);
    //return x0;
    return low;
}

void isolate_root_with_cluster(double& low, double& hgh, Polynomial *p)
{
  dprint(low);
  dprint(hgh);
  // ensure starting conditions
  //
  PointVal pv_low(low,p);
  // double p_low = eval(low);
  if(pv_low.v() < 0.)
  {
    Polynomial p_minus = -(*p);
    p_minus.isolate_root(low,hgh);
    return;
  }

  // initialize algorithm
  //
  // initialize root bracket.
  //
  PointVal pv_hgh(hgh,p);
  RootBracketClusterState rb(pv_low,pv_hgh);
  dprint(rb.low().x());
  dprint(rb.low().v());
  dprint(rb.hgh().x());
  dprint(rb.hgh().v());
  //
  Polynomial Dp = p->get_derivative();

  // requirements of algorithm:
  // - second-order convergence in all cases
  // - high performance in generic cases
  // - handles repeated polynomial roots
  // - converges to a root between low and hgh
  // - residual of answer is less than machine precision

  // conditions maintained by algorithm:
  // * low < hgh
  // * p(low) > EPSILON
  // * p(hgh) < -EPSILON
  // using_hgh = fabs(p(hgh)) > fabs(p(low))

  // state of algorithm is specified by:
  // * low
  // * hgh
  // * using_hgh
  // * num_roots_in_cluster = assumed odd number of roots in cluster

  // iteration of algorithm:
  // use newton iteration on smaller side to look for root
  // - check that guess is between low and hgh
  //   and that new value has smaller residual
  //   else switch to secant method
  // - if newton does not change the sign
  //   increment num_roots_in_cluster until you leave interval,
  //   the sign changes, or the residual increases.
  //  - if the sign does change
  //    decrement num_roots_in_cluster until num_roots_in_cluster=1
  //    or the sign does not change (if residual does not shrink
  //    sufficiently modify num_roots_in_cluster).
  //  - do one iteration of bisection if residual did not
  //    shrink by factor of e.g. .5

  // assumed number of roots in root cluster we are converging towards
  goto start_newton_iteration;

 start_newton_iteration:
  // initialize newton iteration
  // choose smaller endpoint as starting value.
  goto newton_iteration;

 done:
   dprint(rb.num_evals());
   low = rb.low().x();
   hgh = rb.hgh().x();
   //return rb.x();
   return;

 // newton's method begins with a point which is one end of a root bracket.
 // it determines a displacement based on the value of the function
 // and its derivative at the point.  To handle root clusters
 // it assumes that the function is approximated by a shifted
 // and scaled monomial.
 newton_iteration:
  {
    if(rb.done()) goto done;
    dprintf("using newton: x0=%24.16e, p(x0)=%24.16e",rb.pv().x(),rb.pv().v());

    const PointVal start_pv0(rb.pv());
    //
    const double Dpx0 = Dp.eval(rb.pv().x());
    // bail if the derivative is too small
    if(fabs(Dpx0) < 1e-6)
    {
      dprintf("bad derivative: %f",Dpx0);
      goto newton_bailed;
    }
    double displacement = -rb.pv().v()/Dpx0;
   compute_x0:
    double x0 = rb.pv().x() + rb.num_roots_in_cluster()*displacement;
    // if Newton prediction is outside interval
    if(rb.not_in_interval(x0))
    {
      // reset num_roots_in_cluster and try again.
      if(rb.num_roots_in_cluster()!=1)
      {
        rb.reset_num_roots_in_cluster();
        goto compute_x0;
      }
      dprintf("Newton %24.16e outside interval [%24.16e,%24.16e]",
        x0, rb.low().x(), rb.hgh().x());
      dprintf("with values [%24.16e,%24.16e]", rb.low().v(), rb.hgh().v());
      goto newton_bailed;
    }
    // evaluate the polynomial at the new guess
    if(rb.eval(x0,p)) goto done;
    // bail if the new residual is not sufficiently smaller
    if(rb.pv().a() > .5*start_pv0.a()) // exp(-1) < .5 < 1
    {
      dprintf("apx0=%24.16e and start_apx0=%24.16e",rb.pv().a(),start_pv0.a());
      dprintf("so Newton didn't help enough, maybe because Dpx0=%24.16e.",Dpx0);
      goto newton_bailed;
    }
    // did the sign change?
    bool sign_changed = ((start_pv0.v() < 0. && rb.pv().v() > 0.)
                      || (start_pv0.v() > 0. && rb.pv().v() < 0.)) ? true : false;
    // does changing number of roots in root cluster help?
    const PointVal old_pv0 = rb.pv();
    if(sign_changed)
    {
      // keep decrementing presumed number
      // of roots in cluster as long as you can
      // until it does not help the residual
      if(rb.num_roots_in_cluster()>=3)
      {
        //double old_x0 = x0;
       decrement:
        double x0 = old_pv0.x() - 2*displacement;
        // displacement stayed in root bracket, so smaller one will too
        if(rb.eval(x0,p)) goto done;
        if(rb.pv().a() < old_pv0.a())
        {
            // accept decrement
            rb.decrement_num_roots_in_cluster();
            goto decrement;
        }
        else
        {
          // reject decrement
          rb.set_to_using_smaller();
          goto newton_iteration;
        }
      }
      // only one root so just iterate
      goto newton_iteration;
    }
    else // !sign_changed
    {
      // if the residual did not shrink by trigger_ratio then
      // we will try incrementing num_roots_in_cluster
      //
      // compute trigger_ratio
      //
      const double trigger_ratio
        = get_trigger_ratio_for_newton(rb.num_roots_in_cluster());
      // if residual shrunk by more than trigger_ratio we're fine
      if(rb.pv().a() < trigger_ratio*old_pv0.a())
      {
        // evidently not near a root cluster
        // with a wrong assumption about the number
        // of roots, so just iterate
        goto newton_iteration;
      }
      else
      {
        // keep incrementing presumed number of roots
        // in cluster until we are outside the interval
        // or until doing so does not improve
        // the residual by better than multiplying
        // it by a factor of .5
        // (can use any number between 1/e and 1).
        {
         increment:
          x0 += 2*displacement;
          if(rb.not_in_interval(x0))
          {
            // reject increment
            goto newton_iteration;
          }
          const PointVal curr_pv0(rb.pv());
          if(rb.eval(x0,p)) goto done;
          if(rb.pv().a() > .5*curr_pv0.a())
          {
            // reject increment
            //
            // the constraint that we always
            // tighten the root bracket whenever
            // we evaluate the polynomial means
            // we cannot completely reject the evaluation
            // although we can still refrain from
            // incrementing num_roots_in_cluster.
            rb.set_to_using_smaller();
            goto newton_iteration;
          }
          else
          {
            // accept increment
            rb.increment_num_roots_in_cluster();
            goto increment;
          }
        }
      }
      goto newton_iteration;
    }
  }

 newton_bailed:
  {
    // ideally we should try the secant method
    // before reverting to bisection
    //
    // the secant method should handle repeated roots better
    goto secant_bracket_iteration;
  }

 // the secant rule needs to maintain a list of two
 // points that we are "using".
 // the root bracket contains the most recently evaluated point.
 // The secant method is initialized with a root bracket.
 secant_bracket_iteration:
  {
    dprintf("using secant: low=%24.16e, hgh=%24.16e, p_low=%24.16e, p_hgh=%24.16e",
      rb.low().x(),rb.hgh().x(),rb.low().v(),rb.hgh().v());
    // rb.set_to_using_smaller();
    PointVal pv1(rb.pv_old());
    const double start_pva = fmax(pv1.a(),rb.pv().a());

   secant_iteration:
    // double accepted_num_roots_in_cluster = num_roots_in_cluster;
   compute_x0_secant: // for try_num_roots_in_cluster
    const double ratio = rb.pv().v()/pv1.v();
   compute_secant_rule:
      // problem: these numbers can be equal.
      //dprint(rb.pv().v());
      //dprint(pv1.v());
    double x0 = secant_rule(rb.pv().x(),pv1.x(),ratio,
      rb.num_roots_in_cluster());
    // make sure that x0 is in root bracket
    if(x0<rb.low().x() || x0>rb.hgh().x())
    {
      // different signs would force x0 in root bracket
      //assert_gt(ratio, 0.)
      // increasing num roots would only make it worse
      if(rb.num_roots_in_cluster()>1)
      {
        rb.decrement_num_roots_in_cluster();
        goto compute_secant_rule;
      }
      //else
      //{
      //  goto bisect;
      //}
    }
    // evaluate at the candidate
    dprint(x0);
    pv1 = rb.pv(); if(rb.eval(x0,p)) goto done;
    if(rb.pv().a() >= pv1.a())
    {
      dprintf("bisect");
      goto bisect;
    }

    // did residual change signs?

    bool sign_changed = ((pv1.v() < 0. && rb.pv().v() > 0.)
                      || (pv1.v() > 0. && rb.pv().v() < 0.)) ? true : false;
    // if residual changed signs then try decrementing num_roots_in_cluster
    // if residual did not change signs then try incrementing.
    if(!sign_changed)
    {
      // if the residual did not shrink by trigger_ratio then
      // we will try incrementing num_roots_in_cluster
      double secant_trigger = get_trigger_ratio_for_secant(
        rb.num_roots_in_cluster());
      if(pv1.a() > secant_trigger * rb.pv().a())
      {
        // try incrementing num_roots_in_cluster
        const double ratio = rb.pv().v()/pv1.v();
        const double x0 = secant_rule(rb.pv().x(),pv1.x(),ratio,
          rb.num_roots_in_cluster()+2);
        // make sure root is in interval
        if(x0<rb.low().x() || x0>rb.hgh().x())
        {
          dprintf("bisect");
          goto bisect;
        }
        pv1 = rb.pv(); if(rb.eval(x0,p)) goto done;
        if(rb.pv().a() < pv1.a())
        {
          rb.increment_num_roots_in_cluster();
          goto secant_iteration;
        }
      }
    }
    else
    {
      // try decrementing num_roots_in_cluster
      if(rb.num_roots_in_cluster()>=3)
      {
        const double ratio = rb.pv().v()/pv1.v();
        const double x0 = secant_rule(rb.pv().x(),pv1.x(),ratio,
          rb.num_roots_in_cluster()-2);
        // make sure root is in interval
        if(x0<rb.low().x() || x0>rb.hgh().x()) goto bisect;
        pv1 = rb.pv(); if(rb.eval(x0,p)) goto done;
        if(rb.pv().a() < pv1.a())
        {
          rb.decrement_num_roots_in_cluster();
          goto secant_iteration;
        }
      }
    }
   secant_iter_done:
    {
      if(rb.pv().a() > .5*start_pva)
      {
        dprintf("secant did not help enough: apx=%24.16e, old_apx=%24.16e",
          rb.pv().a(), start_pva);
        goto bisect;
      }
      goto secant_bracket_iteration;
    }
  }
  
 bisect:
  {
    dprintf("bisecting: x=%24.16e, px=%24.16e",rb.pv().x(),rb.pv().v());
    const double x0 = (rb.low().x() + rb.hgh().x())*.5;
    if(rb.eval(x0,p)) goto done;
    //goto newton_iteration;
    goto secant_bracket_iteration;
  }
  eprintf("should not get here");
  return;
}

#if 0
  // isolate root using modified secant iteration
  //
  // algorithm:
  //   conditions maintained:
  //   * t_low < t_hgh
  //   * pos(t_low) >= 0
  //   * pos(t_hgh) !>= 0
  //   iteration of algorithm:
  //   * Let t1,t2 be secant and "newton" estimates
  //      for the "worst-offending" indicator
  //      (lie on opposite sides for quadratic)
  //      if I want to sample rather than work with polynomials
  //      explicitly instead of newton I can evaluate at the secant
  //      point and then use it and the endpoint that agrees with it
  //      to get a secant that will lie on the opposite side
  //   * check if conditions are maintained, else revert to modified weighting.
  //   * check if conditions are maintained, else revert to bisection.
  //   * we are done when (t_hgh-t_low) < t_epsilon.
  //
  // idea: move both endpoints each iteration.  goal is to move the big one.
  //   revert to bisection.
  //
  //   possible requirements for worst offendor:
  //   * minimum of secant-rule t values
  //     for indicators that are negative at t_hgh
  //
  double t_hgh=1.;
  bool done=false;
  int miteration=0;
  while(!done)
  {
    // maximum possible number of iterations is
    // (6+4+2) indicator roots per point times numPoints 
    assert_lt(miteration++ < 12*numPoints);
    // identify the worst-offending indicator.
    // The number of roots is finite, so as long as t_hgh
    // moves one root to the left each time this should
    // eventually terminate.
    int mp, mi;
    double p1;
    t_low = 0.; // use regula falsi each time to pick indicator
    double t1 = get_lowest_secant_rule_t(t_low,t_hgh,pos_low,pos_hgh,mp,mi,p1);
    if(!mp)
    {
      // There were no negative points so we are done
      theta = t1;
      break;
    }
    //double p1 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t1);

    // iterate the worst offender to find root.
    // in each iteration of secant method bail and bisect if
    // anything unexpected happens
    //
    p_low = pos_low.get(mp,mi);
    p_hgh = pos_hgh.get(mp,mi);
    int niteration=0;
    while()
    {
      // enforcing interval-halving should ensure termination
      // assert(niteration++ < -log(EPSILON)/log(2));
      assert(niteration++ < 100);
      // use t1 and point whose sign agrees with sign at t1
      // to get a second candidate
      if((t_hgh-t_low)<EPSILON) break;
      const double tL = t_low;
      const double tH = t_hgh;
      assert_lt(tL,t1);
      assert_lt(t1,tH);
      const double pL = p_low;
      const double pH = p_hgh;
      const double tw = t_hgh-t_low;
      // evaluate the indicator at t1
      if(p1<-EPSILON) // t1 will probably be the new tH, so use it:
      {
        t_hgh = t1;
        p_hgh = p1;
        const double dp = pH-p1;
        if(dp>-EPSILON) goto bisect; // expecting negative slope
        const double t2 = (pH*t1-p1*tH)/dp;
        if(t2<=tL.) goto bisect; // expecting t2 to lie in interval
        assert_lt(tL,t2);
        assert_lt(t2,t1);
        double p2 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t2);
        if(p2>EPSILON)
        // p1 and p2 have opposite signs as expected.
        {
          t_low = t2;
          p_low = p2;
          goto get_new_t1;
        }
        else if(p2<-EPSILON)
        // both candidates have the same sign so move the
        // boundary that shares its sign and use an altered
        // weight rule to try to move the other boundary
        {
          t_hgh = t2;
          p_hgh = p2;
          goto get_new_t1_using_half_weight_at_t_low;
        }
        else // we have a root: -EPSILON <= p2 <= EPSILON
        {
          t = t2;
          p = p2;
          goto done;
        }
      }
      else if(p1>EPSILON)
      {
        t_low = t1;
        p_low = p1;
        // t1 will probably be the new t_low, so use it:
        const double dp = p1-pL;
        if(dp>-EPSILON) goto bisect;
        const double t2 = (p1*tL-pL*t1)/dp;
        if(t2>=tH) goto bisect; // expecting t2 to lie in interval
        assert_lt(t2,tH);
        assert_lt(t1,t2);
        double p2 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t2);
        if(p2<-EPSILON)
        // p1 and p2 have opposite signs as expected.
        {
          t_hgh = t2;
          p_hgh = p2;
          goto get_new_t1;
        }
        else if(p2>EPSILON)
        // both candidates have the same sign so move the
        // boundary that shares its sign and use an altered
        // weight rule to try to move the other boundary
        {
          t_low = t2;
          p_low = p2;
          goto get_new_t1_using_half_weight_at_t_hgh;
        }
        else // we have a root: -EPSILON <= p2 <= EPSILON
        {
          t = t2;
          p = p2;
          goto done;
        }
      }
      else // we have a root: -EPSILON <= p1 <= EPSILON
      {
        t = t1;
        p = p1;
        goto done;
      }

      eprintf("should not get here");

      // This should be called whenever the interval does not
      // shrink by a factor of 0.5
      bisect:
      {
        t1 = 0.5*(t_low+t_hgh);
        p1 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t1);
      }

      get_new_t1_using_half_weight_at_t_hgh:
      {
        t1 = (0.5*p_hgh*t_low-p_low*t_hgh)/(0.5*p_hgh-p_low);
        // case p1<0: t1 will become the new t_low (else will bisect anyway)
        // if the width did not shrink enough then bisect
        if((t_hgh-t1)/tw > 0.5) goto bisect;
        p1 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t1);
        // if the sign still did not change then bisect
        if(p1>0.) goto bisect;
        continue;
      }
      //
      get_new_t1_using_half_weight_at_t_low:
      {
        t1 = (p_hgh*t_low-0.5*p_low*t_hgh)/(p_hgh-0.5*p_low);
        // case p1>0: t1 will become the new t_hgh (else will bisect anyway)
        // if the width did not shrink enough then bisect
        if((t1-t_low)/tw > 0.5) goto bisect;
        p1 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t1);
        // if the sign still did not change then bisect
        if(p1<0.) goto bisect;
        continue;
      }

      get_new_t1:
      {
        if((t_hgh-t_low)/tw > 0.5) goto bisect;
        t1 = (p_hgh*t_low-p_low*t_hgh)/(p_hgh-p_low);
        p1 = get_gas10_positivity_indicator_at_point(mp,mi,Q0,dQ,t1);
        continue;
      }

      done:
      {
        break;
      }

}
#endif

void Polynomial::polish_roots(double* roots, int num_roots) const
{
  assert_le(num_roots, order);
  // assume that roots are ordered and represent all real roots
  for(int i=1;i<num_roots;i++){
    assert_le(roots[i-1],roots[i]);
  }
  //
  // we will reject any polishing that causes a root to cross
  // the half-way point between it and neighboring roots.
  //
  double mid_roots[MAX_ORDER+1];
  mid_roots[0] = -DBL_MAX;
  for(int i=1;i<num_roots;i++) mid_roots[i] = 0.5*(roots[i-1]+roots[i]);
  mid_roots[num_roots] = DBL_MAX;

  Polynomial p1 = get_derivative();

  // for each root iterate on newton until the iteration
  // succeeds (the residual becomes small)
  // or fails (the root goes outside its box or the residual becomes larger)
  for(int i=0;i<num_roots;i++)
  {
    //double p0v = pvals[i];
    double x = roots[i];
    double p0v = eval(x);
    int numiter;
    for(numiter=0;numiter<10;numiter++)
    {
      if(fabs(p0v) < 1e-13) break;
      double p1v = p1.eval(x);
      // don't bother trying to polish near excessively clustered roots
      if(fabs(p1v) < 1e-12) // what number should I choose?
      {
        dprint(p1v);
        break;
      }
      x = x - p0v/p1v;
      double p0v_new = eval(x);
      if(fabs(p0v_new)>=fabs(p0v))
      {
        dprintf3("i=%d,numiter=%d,p0v_new=%24.16e,p0v=%24.16e",
           i, numiter, p0v_new, p0v);
        break;
      }
      if(x < mid_roots[i])
      {
        dprintf("i=%d,numiter=%d,x=%24.16e,mid_roots[i]=%24.16e",
           i, numiter, x, mid_roots[i]);
        break;
      }
      if(x > mid_roots[i+1])
      {
        dprintf3("i=%d,numiter=%d,x=%24.16e,mid_roots[i+1]=%24.16e",
           i, numiter, x, mid_roots[i+1]);
        break;
      }
      // accept the polished value.
      dprintf3("i=%d,numiter=%d,p0v=%24.16e,p0v_new=%24.16e",
         i, numiter, p0v, p0v_new);
      p0v = p0v_new;
      roots[i] = x;
    }
  }
}

double Polynomial::polish_root(double x0) const
{
  double x = x0;
  // use one iteration of Newton
  Polynomial p1 = get_derivative();
  double p1v = p1.eval(x);
  if(fabs(p1v) < 1e-14) // what number should I choose?
  {
    // This root is (approximately) a repeated
    // (and hence ill-conditioned) root;
    // Is it repeated only once?
    Polynomial p2 = p1.get_derivative();
    double p2v = p2.eval(x);
    if(fabs(p2v) < 1e-14)
    {
      // the root is (approximately) repeated at least three times
      // approximate it as a cubic?  In any case it is highly
      // ill-conditioned.
      eprintf("implementation incomplete");
    }
    else
    {
      // This root is repeated twice, so find approximate root
      // using the second-order Taylor polynomial.
      eprintf("implementation incomplete");
      // double coef[2];
      // coef[0] = p0v;
      // coef[1] = p1v;
      // coef[2] = p2v;
      // Polynomial p(coef, 2);
      // x -= p.get_smallest_positive_root();
    }
  }
  else
  {
    // find the root of the 1st-order Taylor polynomial
    double p0v = eval(x0);
    x = x - p0v/p1v;
    // should we iterate?
    // for(int i=0; i<num_evals;i++) x = x - eval(x)/p1.eval(x);
  }
  return x;
}

// replace with a better root finder?
// matlab gets a much more accurate answer for
// (x-i)*(x+i)*(x-2)*(x-2) or
// (x-i)*(x+i)*(x-2)
// than this does.
//
// this returns a negative number if no positive root exists
double Polynomial::get_smallest_positive_root() const
{
  // assumes that the leading coefficient is nonzero
  // (e.g. fix_lead_coef() has been applied).
  switch(order)
  {
    case 0:
      invalid_value_error(order);
      return -1.; // NAN
    case 1:
      return -coef[0]/coef[1];
    case 2: // use quadratic formula
    {
      const double& a = coef[0];
      const double bh = coef[1]/2.;
      const double& c = coef[2];
      double discriminant = bh*bh-a*c;
      if(discriminant<0.)
      {
        // a hack to avoid NAN
        if(discriminant > -1e-16)
          discriminant = 0.;
        else
          return -1.; // irreducible quadratic
      }
      const double sq_disc = sqrt(discriminant);
      double t1 = (-bh - sq_disc)/a;
      double t2 = (-bh + sq_disc)/a;
      double tmin,tmax;
      if(t1<t2) {
        tmin = t1;
        tmax = t2;
      }
      else
      {
        tmin = t2;
        tmax = t1;
      }
      if(tmin>=0.) return tmin;
      // if(tmax<0.) return 0./0.;
      return tmax;
    }
    // case 3: use formula for roots of cubic?
    default:
      // reverse order of coefficients
      // (JTrootFinder uses matlab convention)
      double coefficients[MAX_DEGREE+1];
      for(int i=0;i<=order;i++) coefficients[i] = coef[order-i];
      JTrootFinder jtRootFinder;
      double smallest_pos_real_root
        = jtRootFinder.get_smallest_positive_real_root(coefficients, order);
      // should we first polish it?
      return smallest_pos_real_root;
  }
}

class root_bracket
{
  double x_lower;
  double x_upper;
  double f_lower;
  double f_upper;
};

//#include "brent.h"

/* roots/brent.c [EAJ modified]
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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

#define GSL_SUCCESS  0

typedef struct
  {
    double a, b, c, d, e;
    double fa, fb, fc;
  }
brent_state_t;

// in polynomial case we avoid passing function pointer
// to avoid function call overhead.
#define USE_POLYNOMIAL
int brent_init (brent_state_t * state,
#ifdef USE_POLYNOMIAL
  Polynomial *p,
#else
  double (*f)(double x, void * data), void * data,
#endif
  double * root, double x_lower, double x_upper)
{
  *root = 0.5 * (x_lower + x_upper) ;

#ifdef USE_POLYNOMIAL
  const double f_lower=p->eval(x_lower);
  const double f_upper=p->eval(x_upper);
#else
  const double f_lower=f(x_lower,data);
  const double f_upper=f(x_upper,data);
#endif
  
  state->a = x_lower;
  state->fa = f_lower;

  state->b = x_upper;
  state->fb = f_upper;

  state->c = x_upper;
  state->fc = f_upper;

  state->d = x_upper - x_lower ;
  state->e = x_upper - x_lower ;

  if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
    {
      dprint1(x_lower);
      dprint1(x_upper);
      eprintf ("endpoints do not straddle y=0: [%24.16e, %24.16e]\n",
        f_lower, f_upper);
    }

  return 0;

}

/* brent.c -- brent root finding algorithm */
// this returns 1 (not_done) unless the roots
// are very close together or an exact root is found
// (according to machine arithmetic)
int brent_iterate (brent_state_t* state,
#ifdef USE_POLYNOMIAL
  Polynomial *p,
#else
  double (*f)(double x, void * data), void * data,
#endif
  double * root, double * x_lower, double * x_upper, double * f_lower)
{
  double tol, m;

  int ac_equal = 0;

  double a = state->a, b = state->b, c = state->c;
  double fa = state->fa, fb = state->fb, fc = state->fc;
  double d = state->d, e = state->e;
  
  if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
    {
      ac_equal = 1;
      c = a;
      fc = fa;
      d = b - a;
      e = b - a;
    }
  
  if (fabs (fc) < fabs (fb))
    {
      ac_equal = 1;
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
  
  tol = 0.5 * GSL_DBL_EPSILON * fabs (b);
  // tol = 5000 * GSL_DBL_EPSILON * fabs (b);
  // tol = 1e-6 * fabs (b);
  // const double tol_bisection = 1e-10 * fabs (b);
  const double tol_bisection = tol;//10*tol;
  m = 0.5 * (c - b);
  
  if (fb == 0)
    {
      *root = b;
      *x_lower = b;
      *x_upper = b;
      
      return 0;
    }
  
  if (fabs (m) <= tol)
    {
      *root = b;

      if (b < c) 
        {
          *x_lower = b;
          *x_upper = c;
        }
      else
        {
          *x_lower = c;
          *x_upper = b;
        }

      return 0;
    }
  
  if (fabs (e) < tol_bisection || fabs (fa) <= fabs (fb))
    {
      dprintf3("using bisection");
      d = m;            /* use bisection */
      e = m;
    }
  else
    {
      dprintf3("using inverse cubic interpolation");
      double p, q, r;   /* use inverse cubic interpolation */
      double s = fb / fa;
      
      if (ac_equal)
        {
          p = 2 * m * s;
          q = 1 - s;
        }
      else
        {
          q = fa / fc;
          r = fb / fc;
          p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
          q = (q - 1) * (r - 1) * (s - 1);
        }
      
      if (p > 0)
        {
          q = -q;
        }
      else
        {
          p = -p;
        }
      
      if (2 * p < Min (3 * m * q - fabs (tol * q), fabs (e * q)))
        {
          e = d;
          d = p / q;
        }
      else
        {
          dprintf3("interpolation failed, fall back to bisection");
          d = m;
          e = m;
        }
    }
  
  a = b;
  fa = fb;
  
  if (fabs (d) > tol)
    {
      b += d;
    }
  else
    {
      b += (m > 0 ? +tol : -tol);
    }
  
#ifdef USE_POLYNOMIAL
  fb = p->eval(b);
#else
  fb=f(b,data);
#endif

  state->a = a ;
  state->b = b ;
  state->c = c ;
  state->d = d ;
  state->e = e ;
  state->fa = fa ;
  state->fb = fb ;
  state->fc = fc ;
  
  /* Update the best estimate of the root and bounds on each
     iteration */
  
  *root = b;
  
  if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) 
    {
      c = a;
      fc = fa;
    }

  if (b < c)
    {
      *x_lower = b;
      *f_lower = fb;
      *x_upper = c;
    }
  else
    {
      *x_lower = c;
      *f_lower = fc;
      *x_upper = b;
    }

  return 1 ;
}

/* find the root in the interval [lower_bound, upper_bound].  */
void brent_root(
#ifdef USE_POLYNOMIAL
  Polynomial *p,
#else
  double (*f)(double x, void * data), void * data,
#endif
  double& x_lower, double& x_upper)
{
#ifdef USE_POLYNOMIAL
  const int MAX_ITERATIONS = 30;
#else
  const int MAX_ITERATIONS = 30;
#endif
  int not_done = 1;
  int brent_iterations = 0;
  double r;
  //double x_lower = lower_bound;
  //double x_upper = upper_bound;

  brent_state_t state;
#ifdef USE_POLYNOMIAL
  brent_init(&state, p, &r, x_lower, x_upper);
#else
  brent_init(&state, f, data, &r, x_lower, x_upper);
#endif

  double f_lower;
  do 
  {
    brent_iterations++ ;
    dprint3(brent_iterations);

#ifdef USE_POLYNOMIAL
    not_done = brent_iterate(&state, p, &r, &x_lower, &x_upper, &f_lower);
#else
    not_done = brent_iterate(&state, f, data, &r, &x_lower, &x_upper, &f_lower);
#endif
    dprintf3(
      "a =%+24.16e "
      "b =%+24.16e "
      "c =%+24.16e "
      "d =%+24.16e "
      "e =%+24.16e "
      "fa=%+24.16e "
      "fb=%+24.16e "
      "fc=%+24.16e",
      state.a,
      state.b,
      state.c,
      state.d,
      state.e,
      state.fa,
      state.fb,
      state.fc);

    assert_le(x_lower,r);
    assert_le(r,x_upper);

    // double diff_points = fabs(x_upper-x_lower);
    if(not_done)
    {
      //double diff_residual = fabs(state.fb);
      assert_ge(f_lower,0.);
      double diff_residual = f_lower;
      not_done = diff_residual < GSL_DBL_EPSILON ? 0 : 1;
    }
  }
  while (not_done && brent_iterations < MAX_ITERATIONS);

  dprint3(brent_iterations);
  if (brent_iterations == MAX_ITERATIONS)
  {
    // This is taking too long; just use bisection, which
    // requires <46 iterations on interval [0,1] for
    // double precision arithmetic
    //return isolate_root_using_bisection(state.a,state.c,p);
    isolate_root_using_bisection(x_lower,x_upper,p);
    return;
#ifdef USE_POLYNOMIAL
    // taking too long, so check if near a cluster
    //
    // this helps with multiple roots but does not
    // seem to help with clusters; exact multiple
    // roots are presumably rare enough that we
    // can just revert to bisection for that case
    // (which takes about 50 iterations)
    //
    // I was concerned about clusters because
    // in the GEM problem near-gyrotropy generally
    // holds (outside of the diffusion region)
    // where the smaller eigenvalues of the
    // pressure tensor are nearly identical.
    // But in practice when positivity limiters
    // are triggered I find that brent converges
    // in 2 to 4 iterations, so I won't worry about
    // this until I see that I need to. -eaj
    //
    // this code has an error: can get ratio of 1
    // in secant method (change it so it just reverts
    // to bisection)
    //return isolate_root_with_cluster(state.a,state.c,p);
    isolate_root_with_cluster(x_lower,x_upper,p);
    return;
#else
    eprintf("exceeded maximum number of iterations");
#endif
  }
  //return state.b;
  return;
}

// the sign of this function agrees with the sign of the
// polynomial and this is actually a continuous function in a
// topology that includes +inf and -inf (assumes that machine
// properly handles inf); this function cannot have roots with
// multiplicity
double quasi_rational_function(double x, void * data_in)
{
  Polynomial ** data = (Polynomial**) data_in;
  Polynomial *p = data[0];
  Polynomial *dp = data[1];
  double numerator = p->eval(x);
  // 0./0. should be zero (not nan) in this case
  if(numerator==0.) return 0.;
  double denominator = dp->eval(x);
  dprint3(numerator);
  dprint3(denominator);
  double val = numerator/fabs(denominator);
  dprint3(val);
  return val;
}

void Polynomial::isolate_root(double& low, double& hgh)
{
  // brent_root is better in general even though it is slower for
  // roots with multiplicity. multiplicity is rare and for double
  // precision bisection is not a terrible algorithm (maybe
  // would be better to find the roots of the rational function
  // consisting of the polynomial divided by its derivative,
  // since then there are no multiple roots, although eventually
  // the problem becomes ill-conditioned because you are taking
  // the ratio of small numbers.)
  #if 0
  Polynomial dp = get_derivative();
  Polynomial* data[2];
  data[0] = this;
  data[1] = &dp;
  #endif
  //return isolate_root_using_bisection(low,hgh,this);
#ifdef USE_POLYNOMIAL
  brent_root(this, low, hgh);
  return;
#else
  // In practice using a quasi-rational function to
  // avoid rational roots generally doubles the number
  // of polynomial evaluations and even in the case
  // of a root with multiplicity it does not decrease
  // the number of iterations enough to make it faster
  // than bisection.  (The problem is that division
  // of small numbers in finite-precision arithmetic
  // is ill-conditioned.)
  brent_root(quasi_rational_function, data, low, hgh);
  return;
#endif
  // This wastes too much time trying to guess the the number
  // of roots in the cluster. It is faster for roots with
  // multiplicity, but for clusters with nonzero separation when
  // the algorithm gets inside the cluster it has to detect a
  // change in the number of roots which takes time, so the
  // algorithm is not faster in this much more common case.
  // If you have repeated roots you ultimately need to revert to
  // bisection, which is actually not so bad for repeated roots
  // if all you want is a sufficiently small residual (and it's
  // hard to require more than that).
  isolate_root_with_cluster(low,hgh,this);
  return;
}

