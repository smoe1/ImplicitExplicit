#include "dog_math.h"
#include "dogdefs.h"
#include "DogParamsCart2.h"

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
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
  //const int ksize = (method1==2)?1:3;
  int ksize;
  if (method1==2)
    {  ksize = 1;  }
  else
    {  ksize = 3;  }
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
  
  // Moment limiter (highest to lowest)
  switch ( method1 )
    {
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
	      double dwp        = dw_right.get(i,j,m,1);
	      double dwm        =  dw_left.get(i,j,m,1);
	      double wx_now     =  wx_cent.get(i,j,m,1);
	      double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
              
	      wx_cent.set(i,j,m,1, wx_limited );
	      
	      // y-direction: limit linear term
	      dwp        =    dw_up.get(i,j,m,1);
	      dwm        =  dw_down.get(i,j,m,1);
	      double wy_now     =  wy_cent.get(i,j,m,1);
	      double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
	      
	      wy_cent.set(i,j,m,1, wy_limited );
	    }

      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      break;

    case 3:  // 3rd order in space
      
      // for this limiter there is actually no reason to
      // transform the higher-order terms.  This should
      // be changed to accelerate the code.
      //
      // Convert to characteristic variables in x-direction
      ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);
      //
      // Convert to characteristic variables in y-direction
      ConvertQtoW(2,aux,q,dw_up,dw_down,wy_cent,ProjectLeftEig);
      
#undef report_stats
// so that we don't waste computation time if we won't report
#ifdef report_stats
  #define increment(a) a++
#else
  #define increment(a) (void)0
#endif
      const double dx=dogParamsCart2.get_dx();
      const double dy=dogParamsCart2.get_dy();
      // this should be made a configurable parameter
      const double M=0.; //.01; //M=.000050;
      const double Mdx2=M*dx*dx;
      const double Mdy2=M*dy*dy;
#ifdef report_stats
// the accumulation of these variables is not parallelized
      int nx_steep = 0;
      int nxL1 = 0;
      int nxL_mxd = 0;
      int nxL_quad = 0;
      int nxL_1and2 = 0;
      int n = 0;

      int ny_steep = 0;
      int nyL1 = 0;
      int nyL_mxd = 0;
      int nyL_quad = 0;
      int nyL_1and2 = 0;
#endif
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)	
	for (int j=(3-mbc); j<=(my+mbc-2); j++)	      
	  for (int m=1; m<=meqn; m++)
	    {
              increment(n);
              // limit terms in x direction
              bool limited_linear = false;
              if(1)
              {
		// x-direction: limit linear term
                bool limited_linear_x=false;
                if(1)
                {
		  double wx_now     =  wx_cent.get(i,j,m,1);
                  if(fabs(wx_now) > Mdx2) // so the peaks aren't clipped
                  {
                    increment(nx_steep);
		    double dwp        = dw_right.get(i,j,m,1);
		    double dwm        =  dw_left.get(i,j,m,1);
		    double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
                    if(wx_now!=wx_limited)
                    {
                      increment(nxL1);
                      limited_linear_x=true;
                      limited_linear=true;
		      wx_cent.set(i,j,m,1, wx_limited );
                      // zero out higher-order terms
                      //{
		      //  wx_cent.set(i,j,m,2, 0. );
		      //  wx_cent.set(i,j,m,3, 0. );
                      //}
                    }
                  }
                }
              }
	      
              // limit terms in y direction
              if(1)
              {
		// y-direction: limit linear term
                bool limited_linear_y = false;
                if(1)
                {
		  double wy_now     =  wy_cent.get(i,j,m,1);
                  if(fabs(wy_now) > Mdy2) // so the peaks aren't clipped
                  {
                    increment(ny_steep);
		    double dwp        =    dw_up.get(i,j,m,1);
		    double dwm        =  dw_down.get(i,j,m,1);
		    double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
                    if(wy_now!=wy_limited)
                    {
                      increment(nyL1);
                      limited_linear_y = true;
                      limited_linear = true;
		      wy_cent.set(i,j,m,1, wy_limited );
                      // zero out higher-order terms
                      //{
		      //  wy_cent.set(i,j,m,2, 0. );
		      //  wy_cent.set(i,j,m,3, 0. );
                      //}
                    }
                  }
                }

                if(limited_linear)
                {
		  wx_cent.set(i,j,m,2, 0. );
		  wx_cent.set(i,j,m,3, 0. );
		  wy_cent.set(i,j,m,2, 0. );
		  wy_cent.set(i,j,m,3, 0. );
                }
	      }
	    }
#ifdef report_stats
      if(nx_steep > 0)
      {
        printf("n=%4d, nx_steep=%4d, nxL1=%4d"
               "\n",
          n, nx_steep, nxL1);
      }
      if(ny_steep > 0)
      {
        printf("n=%4d, ny_steep=%4d, nyL1=%4d"
               "\n",
          n, ny_steep, nyL1);
      }
#endif
           	    
      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);

      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      break;          
    }
  
}
