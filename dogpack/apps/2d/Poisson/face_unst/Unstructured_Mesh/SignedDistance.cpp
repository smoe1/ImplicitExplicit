#include "meshdefs.h"
#include "dog_math.h"

//  Signed distance function: 
//
//      SignedDistance(x,y) < 0 inside the region
//      SignedDistance(x,y) = 0 on the boundary
//      SignedDistance(x,y) > 0 outside of the region
//
double SignedDistance(point pt)
{
  double x = pt.x;
  double y = pt.y;

  double out1 =  (x+1.0);    
  double out2 = -(x-1.0);
  double out3 =  (y+1.0);  
  double out4 = -(y-1.0);    
  double lft_eye = sqrt(pow(x+0.5,2) + pow(y-0.5,2)) - 0.25;
  double rgt_eye = sqrt(pow(x-0.5,2) + pow(y-0.5,2)) - 0.25;
  double mouth = sqrt(pow(x,2) + 10.0*pow(y+0.5,2)) - 0.75;
  double nose = sqrt(3.0*pow(x,2) + pow(y-0.1,2)) - 0.15;
    
  double dist = -Min(Min(Min(Min(out1,out2),
			     Min(out3,out4)),
			 Min(lft_eye,rgt_eye)),
		     Min(mouth,nose));
  
  return dist;
}
