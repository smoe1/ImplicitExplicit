function alpha = GetAlpha( rho, p, q, r )

    CC=[13.0*p^3, -6.0*p^3, p*(-12.0*p^2 + 5.0*r*rho), - 5.0*rho* ...
        q^2];
    
    MR = roots(CC);
    
    for k=1:3
        if abs(imag(MR(k)))<1.0e-13 && real(MR(k))>=3.0/13.0 && real(MR(k))<=1.0
            alpha = real(MR(k));
        end
    end
    
end