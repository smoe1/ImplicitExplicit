
def GetAlpha( rho, p, q, r ):

    import numpy as np
    from math import sqrt
    from math import fabs
    from math import pi
    from math import acos
    from math import cos

    if fabs(q) <= 1.0e-12:
      if r<3.0*pow(p,2)/rho:
         return sqrt((3.0*pow(p,2)-rho*r)/(2*pow(p,2)))
      else:
         return 0.0
      exit()
    
    else:
      p2 = p*p;
      p3 = p*p2;
      q2 = q*q;
      Q  = 1.0/6.0*(3.0-rho*r/p2)
      R  = rho*q2/(4.0*p3)
      R2 = R*R
      Q3 = pow(Q,3)
      D  = R2-Q3            
      
      if D >= 0.0:
        
	  SD = sqrt(D)
          Tstar = R - SD
          if Tstar < 0.0:
           
              T = -pow(fabs(Tstar),1.0/3.0)
           
          else: 
            
              T = pow(fabs(Tstar),1.0/3.0)
            
          return pow(R+SD,1.0/3.0)+T
        
      else:
        
	  SQ = sqrt(Q)
          theta = 1.0/3.0 * acos(R/sqrt(Q3))
          a1 = 2.0*SQ*cos(theta)
          if  a1>=0.0 and a1<=1.0:
            
              return a1            
            
          else:
            
              a2 = 2.0*SQ*cos(theta+2.0/3.0*pi);
              if a2>=0.0 and a2<=1.0:
                
                  return a2                 
                
              else:
                
                  a3 = 2.0*SQ*cos(theta+4.0/3.0*pi)
                  return a3                  
                
            
        
    
