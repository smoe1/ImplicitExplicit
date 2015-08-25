function [ dfout ] = DZeqns( zin )

    z1 = zin(1);
    z2 = zin(2);
    z3 = zin(3);
    z4 = zin(4);
    z5 = zin(5);
    
    dfout = zeros(5,5);
    
    dfout(1,1) = z2*z3*z4*z5;
    dfout(1,2) = z1*z3*z4*z5;
    dfout(1,3) = z1*z2*z4*z5;
    dfout(1,4) = z1*z2*z3*z5;
    dfout(1,5) = z1*z2*z3*z4;
    
    dfout(2,1) = -z4*z3*z2 - z2*z4*z5 - z2*z3*z5 - z3*z4*z5;
    dfout(2,2) = -z1*z3*z4 - z4*z5*z1 - z1*z3*z5 - z3*z4*z5;
    dfout(2,3) = -z4*z2*z1 - z1*z2*z5 - z2*z4*z5 - z4*z5*z1;
    dfout(2,4) = -z3*z2*z1 - z1*z2*z5 - z1*z3*z5 - z2*z3*z5;
    dfout(2,5) = -z4*z2*z1 - z3*z2*z1 - z4*z3*z2 - z1*z3*z4;
                          
    dfout(3,1) = z2*z4 + z3*z2 + z2*z5 + z3*z4 + z3*z5 + z4*z5;
    dfout(3,2) = z4*z1 + z3*z1 + z5*z1 + z3*z4 + z3*z5 + z4*z5;
    dfout(3,3) = z2*z1 + z4*z1 + z5*z1 + z2*z4 + z2*z5 + z4*z5;
    dfout(3,4) = z5*z1 + z2*z1 + z3*z1 + z2*z5 + z3*z5 + z3*z2;
    dfout(3,5) = z2*z1 + z3*z1 + z4*z1 + z2*z4 + z3*z2 + z3*z4;
    
    dfout(4,1) = -z2 - z4 - z3 - z5;
    dfout(4,2) = -z1 - z4 - z3 - z5;
    dfout(4,3) = -z2 - z1 - z4 - z5;
    dfout(4,4) = -z1 - z2 - z3 - z5;
    dfout(4,5) = -z2 - z1 - z4 - z3;
    
    dfout(5,1) = 1;
    dfout(5,2) = 1;
    dfout(5,3) = 1;
    dfout(5,4) = 1;
    dfout(5,5) = 1;

end