%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%                   meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc): grid points (cell centers), size = (mx,my)
%%       (xl,yl): grid points (lower left cell corners), size = (mx+1,my+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clear qex;
  qex = zeros(mx,my);
  t   = time;
  for j=1:my
    for i=1:mx
      x = xc(i,j);   v = yc(i,j);
%     qex(i,j) = ( 2 - cos( 2*(x-pi*t) ) )*exp(-4*v^2); 
      qex(i,j) = ( 2 - cos( 2*(x-pi*t) ) )*exp(-4*(v-0.25)^2); 
    end
  end

%% disp([' time = ', num2str(time) ] );

  err = reshape(abs(qex-qsoln(:,:,1)),mx*my,1);
  err_scale = reshape(abs(qex),mx*my,1);
  
  err_rel = norm(err,2)/norm(err_scale,2);
  
  disp([num2str(err_rel,'%1.3e')]);
