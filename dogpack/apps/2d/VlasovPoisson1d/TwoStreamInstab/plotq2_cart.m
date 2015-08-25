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

  figure(1);
  clf;

  %%% linux machines don't choke with contourf like they do with surf %%%
% surf(xl,yl,qaug(:,:,m));
  pcolor(xl,yl,qaug(:,:,m));

  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis([-2*pi 2*pi -2*pi 2*pi]);
  set(gca,'fontsize',16);
  t1 = title(['f(t,x,v) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  shading flat;
  colorbar;
  caxis auto;
% caxis([0,0.37]);
  set(gca,'fontsize',16);
  set(gca,'xtick',-6:2:6);
  set(gca,'ytick',-6:2:6);

  disp([['  minimum value of f is = ', num2str(min(min(qsoln)), '%2.8e')]])

  figure(2)
  clf;
  if mod(mx,2) == 0 % even number of cells printed!
    xcenter_index = mx / 2;
  else
    xcenter_index = (mx+1) / 2;
  end

  vel = yc(xcenter_index,:);
  qslice = qsoln(xcenter_index,:)';
  pp2=plot(vel, zeros(size(vel)), 'b-');  % show y = 0 on the grid
  set(pp2,'linewidth', 2.5);
  hold on;
  pp1=plot(vel, qslice,'r-');
  set(pp1,'linewidth', 2.5);
  axis([yl(1,1) yl(1,my) -0.05 0.45]);
  hold on;
  
  set(gca,'fontsize',16);
  %t1 = title(['f(x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v,t) at t = ',num2str(time),'     [DoGPack]']); 
  t1 = title(['f(t, x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v) at t = ',num2str(time)]); 
  set(t1,'fontsize',16);

  figure(1);

  %%%% Print the pretty pictures!
  picname1 = strcat(outputdir, '.eps');
% picname2 = strcat(outputdir, '_slice', '.eps');
% picname3 = strcat(outputdir, '_slice2', '.eps');
  print(1, '-depsc2', picname1);
% print(2, '-depsc2', picname2);
% print(3, '-depsc2', picname3);
