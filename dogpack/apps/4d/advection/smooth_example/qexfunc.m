function q = qexfunc(time,x,y,z,params)

    x0             = params(1);
    y0             = params(2);
    z0             = params(3);
    width          = params(4);
    half_pi_owidth = params(5);
    u              = params(6);
    v              = params(7);
    w              = params(8);
    xlow           = params(9);
    xhigh          = params(10);
    ylow           = params(11);
    yhigh          = params(12);
    zlow           = params(13);
    zhigh          = params(14);
    
    xx = ((x-u*time)-xlow)/(xhigh-xlow);
    yy = ((y-v*time)-ylow)/(yhigh-ylow);
    zz = ((z-w*time)-zlow)/(zhigh-zlow);
    
    xs = xlow + (xx-floor(xx))*(xhigh-xlow);
    ys = ylow + (yy-floor(yy))*(yhigh-ylow);
    zs = zlow + (zz-floor(zz))*(zhigh-zlow);
    
    r = sqrt((xs-x0)^2+(ys-y0)^2+(zs-z0)^2);
    if (r<width)
        q = (cos(half_pi_owidth*r))^6;
    else
        q = 0;
    end

end