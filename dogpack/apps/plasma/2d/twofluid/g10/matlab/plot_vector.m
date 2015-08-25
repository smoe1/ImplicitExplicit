% flip = 1 for vector, -1 for pseudo-vector
function plot_vector(xc_in,yc_in,vecx_in,vecy_in,flip,scale);

    if(nargin<6)
       scale=.4;
    end
    % set skip so that
    xnum = size(xc_in,1);
    xnum_desired = 20;
    xskip_f = xnum/xnum_desired;
    xskip = ceil(xskip_f);
    xstart = ceil(xskip_f/2.);

    xc=xc_in(xstart:xskip:end,xstart:xskip:end);
    yc=yc_in(xstart:xskip:end,xstart:xskip:end);
    vecx=vecx_in(xstart:xskip:end,xstart:xskip:end);
    vecy=vecy_in(xstart:xskip:end,xstart:xskip:end);

    if(scale~=0)
      % rescale the magnitudes by rms magnitude for display
      vec2 = (vecx.^2 + vecy.^2);
      vecmag = sqrt(vec2);
      rms_mag = sqrt(mean(mean(vec2)));

      vecx = vecx/rms_mag;
      vecy = vecy/rms_mag;
      vecmag = vecmag/rms_mag;
      vec_scale=1./(1.+vecmag);
      vecx = vecx.*vec_scale;
      vecy = vecy.*vec_scale;
    end

    if(flip==1)
      vx = [-vecx(end:-1:1,end:-1:1), vecx(:,end:-1:1); ...
            -vecx(end:-1:1,:),vecx];
      vy = [-vecy(end:-1:1,end:-1:1),-vecy(:,end:-1:1); ...
             vecy(end:-1:1,:),vecy];
    else
      vx = [-vecx(end:-1:1,end:-1:1),-vecx(:,end:-1:1); ...
             vecx(end:-1:1,:),vecx];
      vy = [-vecy(end:-1:1,end:-1:1), vecy(:,end:-1:1); ...
            -vecy(end:-1:1,:),vecy];
    end
    cx = [-xc(end:-1:1,end:-1:1), xc(:,end:-1:1); ...
          -xc(end:-1:1,:),xc];
    cy = [-yc(end:-1:1,end:-1:1),-yc(:,end:-1:1); ...
           yc(end:-1:1,:),yc];
          
    quiver(cx,cy,vx,vy,scale);

    %quiver(xc,yc,vecx,vecy,scale);
    %if(flip==-1)
    %  quiver(xc,-yc,-vecx,vecy,scale);
    %  quiver(-xc,yc,vecx,-vecy,scale);
    %else
    %  quiver(xc,-yc,vecx,-vecy,scale);
    %  quiver(-xc,yc,-vecx,vecy,scale);
    %end
    %quiver(-xc,-yc,-vecx,-vecy,scale);
end

