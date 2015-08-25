function alpha = GetAlpha( rho, p, q, r )

    if abs(q) <= 1.0e-12
      if r<3.0*p^2/rho
         alpha = sqrt((3.0*p^2-rho*r)/(2*p^2));
      else
         alpha = 0.0;
      end
    
    else
      p2 = p*p;
      p3 = p*p2;
      q2 = q*q;
      Q  = 1.0/6.0*(3.0-rho*r/p2);
      R  = rho*q2/(4.0*p3);
      R2 = R*R;
      Q3 = Q^3;
      D  = R2-Q3;        
      
      if D >= 0.0
        
	  SD = sqrt(D);
          Tstar = R - SD;
          if Tstar < 0.0
           
              T = -abs(Tstar)^(1.0/3.0);
           
          else
            
              T =  abs(Tstar)^(1.0/3.0);
              
          end
            
          alpha = (R+SD)^(1.0/3.0)+T;
        
      else
        
	  SQ = sqrt(Q);
          theta = 1.0/3.0 * acos(R/sqrt(Q3));
          a1 = 2.0*SQ*cos(theta);
          if  a1>=0.0 && a1<=1.0
            
              alpha = a1;
            
          else
            
              a2 = 2.0*SQ*cos(theta+2.0/3.0*pi);
              if a2>=0.0 && a2<=1.0
                
                  alpha = a2;
                
              else
                
                  alpha = 2.0*SQ*cos(theta+4.0/3.0*pi);
                
              end
              
          end
          
      end
      
    end
        
end
