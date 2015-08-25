function [ fout ] = Zeqns( zin, Uin )

    z1 = zin(1);
    z2 = zin(2);
    z3 = zin(3);
    z4 = zin(4);
    z5 = zin(5);
    
    rho = Uin(1);
    p   = Uin(2);
    q   = Uin(3);
    alpha = Uin(4);
    
    rhs  = zeros(5,1);
    fout = zeros(5,1);
    
    rhs(1) = (2*(-1+alpha))*q*p*(2*alpha-3)/(rho^2*alpha);
    rhs(2) = -(-3*rho*q^2-24*p^3*alpha^3+3*rho*alpha*q^2+10*p^3*alpha^4+15*p^3*alpha^2)/(rho^2*p*alpha^2);
    rhs(3) = 2*q*(5*alpha-6)/(rho*alpha);
    rhs(4) = -(rho*q^2+8*p^3*alpha^3-10*p^3*alpha^2)/(rho*p^2*alpha^2);
    rhs(5) = 2*q/(p*alpha);
    
    fout(1) = z1*z2*z3*z4*z5 - rhs(1);
    fout(2) = -(z1*z2*z3*z4+z1*z2*z4*z5+z1*z2*z3*z5+z1*z3*z4*z5+z2*z3*z4*z5) - rhs(2);
    fout(3) = z4*z2*z1+z3*z2*z1+z1*z2*z5+z1*z3*z4+z1*z3*z5+z4*z5*z1+z4*z3*z2+z2*z4*z5+z2*z3*z5+z3*z4*z5 - rhs(3);
    fout(4) = -(z2*z1+z4*z1+z3*z1+z5*z1+z2*z4+z3*z2+z2*z5+z3*z4+z3*z5+z4*z5) - rhs(4);
    fout(5) = z1+z2+z4+z3+z5 - rhs(5);

end

