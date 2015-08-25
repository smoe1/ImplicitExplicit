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
  double xin = pt.x;
  double yin = pt.y;
  double rad = sqrt(pow(xin-0.5,2)+pow(yin-0.5,2));
   
  double dist = rad - 0.5e0;
    
  return dist;
}
