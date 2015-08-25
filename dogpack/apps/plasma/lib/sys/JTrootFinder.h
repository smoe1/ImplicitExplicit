
// ==== begin modification of http://www.crbond.com/download/misc/rpoly.cpp ====
// eaj: I wrapped this up in a class to make it thread-safe
// eaj: could accelerate by making some methods inline
#define MAX_DEGREE 20
class JTrootFinder
{
 public:
  JTrootFinder():stopCode(NO_STOP){}
  int rpoly(double *op, int degree, double *zeror, double *zeroi, int info[] );
  double get_smallest_positive_real_root(double *op, int degree);
  double polish_real_root(double seed);
 private:
  // called twice
  void quad(double a,double b1,double c,double *sr,double *si,
          double *lr,double *li);
  // called once
  void fxshfr(int l2, int *nz);
  void quadit(double *uu,double *vv,int *nz);
  void realit(double sss, int *nz, int *iflag);
  // called many times
  void calcsc(int *type);
  void nextk(int *type);
  // called twice
  void newest(int type,double *uu,double *vv);
  // called 5 times, short
  void quadsd(int n,double *u,double *v,double *p,double *q,
          double *a,double *b);
 private:
  // double *p,*qp,*k,*qk,*svk;  // eaj: changed to static memory allocation
  double sr,si,u,v,a,b,c,d,a1,a2;
  double a3,a6,a7,e,f,g,h,szr,szi,lzr,lzi;
  double eta,are,mre;
  int n,nn,nmi,zerok;
  /*static*/ int itercnt;
  // eaj insert: for performance avoid dynamic memory allocation
  //
  double p   [MAX_DEGREE+1];
  double qp  [MAX_DEGREE+1];
  double k   [MAX_DEGREE+1];
  double qk  [MAX_DEGREE+1];
  double svk [MAX_DEGREE+1];
  //
  enum StopCode {
    NO_STOP=0,
    POSITIVE_REAL_ROOT=1
    //,REAL_ROOT=2 // should support this
  };
  int stopCode; // eaj stop condition
  double desired_root;
};
// ==== end modification of http://www.crbond.com/download/misc/rpoly.cpp ====

