function [z,NumIts,mfound] = NewtonSys(rho,u,p,q,r)

    TOL = 1.0e-13;
    MaxIts = 1000;
    NumIts = 0;
    mstop = 0;
    mfound = 0;
    
    alpha = GetAlpha(rho,p,q,r);
    
    if (alpha<1.0e-5)
        qt = 0;
    else
        qt = q/alpha;
    end
    
    uvec = u*ones(5,1);
    Uin(1) = rho;
    Uin(2) = p;
    Uin(3) = q;
    Uin(4) = alpha;
    
    MGL = 2.856970013872806;
    MGS = 1.355626179974266;
  
    z1guess = 0.4*qt/p - MGL*sqrt((1.0-alpha)*p/rho);
    z2guess = 0.4*qt/p - MGS*sqrt((1.0-alpha)*p/rho);
    z3guess = 0.4*qt/p;
    z4guess = 0.4*qt/p + MGS*sqrt((1.0-alpha)*p/rho);
    z5guess = 0.4*qt/p + MGL*sqrt((1.0-alpha)*p/rho);
    
    z = [z1guess;z2guess;z3guess;z4guess;z5guess];

    Fval = Zeqns(z,Uin);
    if norm(Fval,2)<=TOL
        mstop=1;
    end
    
    while(mstop==0)
        
        NumIts = NumIts + 1;
       
        DFval = DZeqns(z);
        
        zcorr = DFval\Fval;
        z = z - zcorr;
        
        Fval = Zeqns(z,Uin);
        
        if (norm(zcorr,2)<=TOL && norm(Fval,2)<=TOL)
            mstop = 1;
            mfound = 1;
            z = uvec + z;
        elseif (NumIts>=MaxIts)
            mstop = 1;
            disp(['   rho = ',num2str(rho,'%10.5e')]);
            disp(['     p = ',num2str(p,'%10.5e')]);
            disp(['     q = ',num2str(q,'%10.5e')]);
            disp([' alpha = ',num2str(alpha,'%10.5e')]);
            disp(' ');
            disp([' norm(zcorr,2) = ',num2str(norm(zcorr,2),'%10.5e')]);
            disp(['  norm(Fval,2) = ',num2str(norm(Fval,2),'%10.5e')]);
            disp(' ');
            %error(' Not converged ...');
        end
        
    end
            
    
end