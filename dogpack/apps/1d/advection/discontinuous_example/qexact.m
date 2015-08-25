function q = qexact( x )
% Exact solution to the problem

    % parameters used for exact solution:
    a = 0.5;
    z = -0.7;
    delta = 0.005;
    alpha = 10.0;
    beta  = log(2.0)/(36.*delta*delta);

    q = 0.*x;

    % first bell:
    I = find( -0.8 <= x & x <= -0.6 );
    q(I) = ( G(x(I),z-delta,beta) + G(x(I),z+delta,beta) + 4.0*G(x(I),z,beta) ) / 6.0;

    % second bell:
    I    = find( -0.4 <= x & x <= -0.2 );
    q(I) = 1.0 + 0.*x(I);

    % third bell:
    I    = find( 0.0 <= x & x <= 0.2 );
    q(I) = 1.0-abs(10*(x(I)-0.1));

    % fourth bell:
    I    = find( 0.4 <= x & x <= 0.6 );
    q(I) = ( F(x(I),a-delta,alpha) + F(x(I),a+delta,alpha) + 4.0*F(x(I),a,alpha) ) / 6.0;

end

function g = G( x, xc, beta )
    g = exp( -beta*(x-xc).^2 );
end

function f = F( x, xc, alpha )
    f = sqrt( max( 1.0-alpha^2*(x-xc).^2, 0. ) );
end

