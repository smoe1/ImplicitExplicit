#include "dog_math.h"
#include "debug.h"
#include "constants.h"
#include "tensors.h"
//#include "DogParams.h"
//#include "DogParamsCart2.h"
//#include "GEMparams.h"
//#include "GEM.h"
//#include "Legendre2d.h"
//#include <cart/Limiters.h>
//#include "Limiters.h"

// moved here from apps/plasma/lib/2d/ApplyLimiter_divfree.cpp
// (which I copied and modified from lib/2d/cart; I changed
// the coefficients to limiter more aggressively.)
//
void ApplyLimiterKrivodonova_modified(dTensorBC4& aux, dTensorBC4& q,
                  void (*ProjectRightEig)(int,const dTensor1&,
                                          const dTensor1&,const dTensor2&,
                                          dTensor2&),
                  void (*ProjectLeftEig)(int,const dTensor1&,
                                         const dTensor1&,const dTensor2&,
                                         dTensor2&))
{
  // -------------------------------------------------------------  
  // Limiter based on a modification of the following paper:
  //      L. Krivodonova. "Limiters for high-order discontinuous
  //      Galerkin methods." J. Comp. Phys., Vol. 226, pg 879-896.
  // ------------------------------------------------------------- 
  const int   mx = q.getsize(1);
  const int   my = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(2);
  const int method1 = int((sqrt(1+8*kmax)-1)/2);
  const int ksize = (method1==2) ? 1 : 3;
  dTensorBC4  wx_cent(mx,my,meqn,ksize,mbc);
  dTensorBC4  wy_cent(mx,my,meqn,ksize,mbc);
  dTensorBC4 dw_right(mx,my,meqn,ksize,mbc);
  dTensorBC4  dw_left(mx,my,meqn,ksize,mbc);
  dTensorBC4    dw_up(mx,my,meqn,ksize,mbc);
  dTensorBC4  dw_down(mx,my,meqn,ksize,mbc);
  void ConvertQtoW(int ixy, const dTensorBC4& aux, const dTensorBC4& qold, 
                   dTensorBC4& dwp, dTensorBC4& dwm, dTensorBC4& w_cent,
                   void (*ProjectLeftEig)(int,const dTensor1&,
                                          const dTensor1&,const dTensor2&,
                                          dTensor2&));
  void ConvertWtoQ(int ixy, const dTensorBC4& aux, const dTensorBC4& qold,
                   const dTensorBC4& win, dTensorBC4& qout,
                   void (*ProjectRightEig)(int,const dTensor1&,
                                           const dTensor1&,const dTensor2&,
                                           dTensor2&));
  
  // ----------------------------------------------------------
  // Key: storage of characteristic variables
  //       
  //      wx_cent(i,j,me,1)  = Lx * q(i,j,me,2)
  //      wx_cent(i,j,me,2)  = Lx * q(i,j,me,4)
  //      wx_cent(i,j,me,3)  = Lx * q(i,j,me,5)
  //
  //      dw_right(i,j,me,1) = Lx * (q(i+1,j,me,1)-q(i,j,me,1))
  //      dw_right(i,j,me,2) = Lx * (q(i+1,j,me,2)-q(i,j,me,2))
  //      dw_right(i,j,me,3) = Lx * (q(i,j+1,me,2)-q(i,j,me,2))
  //
  //      dw_left(i,j,me,1)  = Lx * (q(i,j,me,1)-q(i-1,j,me,1))
  //      dw_left(i,j,me,2)  = Lx * (q(i,j,me,2)-q(i-1,j,me,2))
  //      dw_left(i,j,me,3)  = Lx * (q(i,j,me,2)-q(i,j-1,me,2))
  //
  //      wy_cent(i,j,me,1)  = Ly * q(i,j,me,3)
  //      wy_cent(i,j,me,2)  = Ly * q(i,j,me,4)
  //      wy_cent(i,j,me,3)  = Ly * q(i,j,me,6)
  //
  //      dw_up(i,j,me,1)    = Ly * (q(i,j+1,me,1)-q(i,j,me,1))
  //      dw_up(i,j,me,2)    = Ly * (q(i+1,j,me,3)-q(i,j,me,3))
  //      dw_up(i,j,me,3)    = Ly * (q(i,j+1,me,3)-q(i,j,me,3))
  //
  //      dw_down(i,j,me,1)  = Ly * (q(i,j,me,1)-q(i,j-1,me,1))
  //      dw_down(i,j,me,2)  = Ly * (q(i,j,me,3)-q(i-1,j,me,3))
  //      dw_down(i,j,me,3)  = Ly * (q(i,j,me,3)-q(i,j-1,me,3))
  //
  // ---------------------------------------------------------- 
  
  // apply Krivodonova limiter to suppress oscillations
  //
  // Moment limiter (highest to lowest)
  switch ( method1 )
    {
    default: unsupported_value_error(method1);
    case 2:  // 2nd order in space   
      
      // Convert to characteristic variables in x-direction
      ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);
      
      // Convert to characteristic variables in y-direction
      ConvertQtoW(2,aux,q,dw_up,dw_down,wy_cent,ProjectLeftEig);
      
      // Limit in both the x-direction and the y-direction
      #pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)          
        for (int j=(3-mbc); j<=(my+mbc-2); j++)       
          for (int m=1; m<=meqn; m++)
            {
              // x-direction: limit linear term
              {
                const double        dwp = dw_right.get(i,j,m,1);
                const double        dwm =  dw_left.get(i,j,m,1);
                const double     wx_now =  wx_cent.get(i,j,m,1);
                const double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
                
                wx_cent.set(i,j,m,1, wx_limited );
              }
              
              // y-direction: limit linear term
              {
                const double dwp        =    dw_up.get(i,j,m,1);
                const double dwm        =  dw_down.get(i,j,m,1);
                const double wy_now     =  wy_cent.get(i,j,m,1);
                const double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
              
                wy_cent.set(i,j,m,1, wy_limited );
              }
            }
      
      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      break;
      
    case 3:  // 3rd order in space
      
      // Convert to characteristic variables in x-direction
      ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);
      
      // Convert to characteristic variables in y-direction
      ConvertQtoW(2,aux,q,dw_up,dw_down,wy_cent,ProjectLeftEig);

      // Limit in the x-direction
      #pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
        for (int j=(3-mbc); j<=(my+mbc-2); j++)   
          for (int m=1; m<=meqn; m++)
            {
              bool limited_x2=false;
              // x-direction: limit quadratic term
              {
                const double dwp        = dw_right.get(i,j,m,2);
                const double dwm        =  dw_left.get(i,j,m,2);
                const double wx_now     =  wx_cent.get(i,j,m,3);
//              const double wx_limited = minmod(wx_now, 0.5*(sq3/sq5)*dwp, 
//                                          0.5*(sq3/sq5)*dwm);
                const double wx_limited = minmod(wx_now, (0.5/sq3)*dwp, 
                                            (0.5/sq3)*dwm);
                if(wx_limited!=wx_now)
                {
                  limited_x2=true;
                  wx_cent.set(i,j,m,3, wx_limited );
                }
              }
              
              //if(limited_x2)
              // x-direction: limit mixed term
              { 
                const double dwp        = dw_right.get(i,j,m,3);
                const double dwm        =  dw_left.get(i,j,m,3);
                const double wx_now     =  wx_cent.get(i,j,m,2);
                const double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
                
                if(wx_limited!=wx_now)
                {
                  limited_x2=true;
                  wx_cent.set(i,j,m,2, wx_limited );
                }
              } 
                  
              if(limited_x2)
              // x-direction: limit linear term
              { 
                const double dwp        = dw_right.get(i,j,m,1);
                const double dwm        =  dw_left.get(i,j,m,1);
                const double wx_now     =  wx_cent.get(i,j,m,1);
                const double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
                
                if(wx_limited!=wx_now)
                {
                  wx_cent.set(i,j,m,1, wx_limited );
                }
              }
                      
              bool limited_y2=false;
              // y-direction: limit quadratic term
              {
                const double dwp    =    dw_up.get(i,j,m,3);
                const double dwm    =  dw_down.get(i,j,m,3);
                const double wy_now =  wy_cent.get(i,j,m,3);
//              const double wy_limited = minmod(wy_now, 0.5*(sq3/sq5)*dwp, 
//                                         0.5*(sq3/sq5)*dwm);
                const double wy_limited = minmod(wy_now, (0.5/sq3)*dwp, 
                                           (0.5/sq3)*dwm);
                if(wy_limited!=wy_now)
                {
                  limited_y2=true;
                  wy_cent.set(i,j,m,3, wy_limited );
                }
              }
              
              //if(limited_y2)
              // y-direction: limit mixed term
              {
                const double dwp        =    dw_up.get(i,j,m,2);
                const double dwm        =  dw_down.get(i,j,m,2);
                const double wy_now     =  wy_cent.get(i,j,m,2);
                const double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
                
                if(wy_limited!=wy_now)
                {
                  limited_y2=true;
                  wy_cent.set(i,j,m,2, wy_limited);
                }
              }
                
              if(limited_y2)
              // y-direction: limit linear term
              {
                const double dwp        =    dw_up.get(i,j,m,1);
                const double dwm        =  dw_down.get(i,j,m,1);
                const double wy_now     =  wy_cent.get(i,j,m,1);
                const double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
                
                if(wy_limited!=wy_now)
                {
                  wy_cent.set(i,j,m,1, wy_limited );
                }
              }
            }

      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      break;          
    }
}
