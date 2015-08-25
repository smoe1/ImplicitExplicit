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
  t = time;
  for j=1:my
    for i=1:mx
      x = xc(i,j);   v = yc(i,j);
%     qex(i,j) = ( 2 - cos( 2*(x-pi*t) ) )*exp(-4*v^2); 
      qex(i,j) = ( 2 - cos( 2*(x-pi*t) ) )*exp(-4*(v-0.25)^2); 
    end
  end

  figure(1);
  clf;
  surf(xl,yl,qaug(:,:,m));
  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis auto;
  axis([xl(1,1) xl(mx,1) yl(1,1) yl(1,my)]);
  set(gca,'fontsize',16);
  t1 = title(['f(x,v,t) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  caxis auto;
  shading flat;
  colorbar;

  err = norm(qsoln-qex,2) / norm(qex,2);

  disp([['  err = ', num2str(min(min(err)), '%2.8e'), '   mx*points_per_cell = ', int2str(mx)]])

  %% Plot err = (qsoln - qex) %%
  figure(2)
  clf;
  surf(xc,yc,qsoln-qex);
  shading flat
  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis auto;
  axis([xl(1,1) xl(mx,1) yl(1,1) yl(1,my)]);
  set(gca,'fontsize',16);
  t1 = title(['qsoln-qex at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  caxis auto;
  shading flat;
  colorbar;


