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
  pcolor(xl,yl,qaug(:,:,m));
  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis auto;
  axis([xl(1,1) xl(mx,1) yl(1,1) yl(1,my)]);
  set(gca,'fontsize',16);
% t1 = title(['f(x,v,t) at t = ',num2str(time),'     [DoGPack]']); 
  t1 = title(['f(t,x,v) at t = ',num2str(time)]); 
  set(t1,'fontsize',16);

  shading flat;
  colorbar;
  caxis auto;
  caxis([0,0.6]);
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
  plot(vel, qslice,'r-');
  axis([yl(1,1) yl(1,my) -0.05 0.65]);
  hold on;
  plot(vel, zeros(size(vel)), 'b-');  % show y = 0 on the grid
  set(gca,'fontsize',16);
  t1 = title(['f(t, x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v) at t = ',num2str(time)]); 
% t1 = title(['f(x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v,t) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);

% figure(3)
% clf;
% E = aux(:,1,2);
% plot(xc, E,'k');
% axis([xl(1,1) xl(mx,1) -0.80 0.80]);
% hold on;
% plot(xc(:,1), zeros(size(xc)), 'b-');  % show y = 0 on the grid
% set(gca,'fontsize',16);
% t1 = title(['E(t,x) at t = ',num2str(time)]); 
% set(t1,'fontsize',16);


% figure(4);
% clf;

% vel = yc(xcenter_index,:);
% I = find( vel > 4.8 );
% vel = vel(I);

% qslice = qsoln(xcenter_index,:)';

% qslice = qslice(I);
% plot(vel, qslice,'r-', 'linewidth', 3);

% hold on;
% plot(vel, zeros(size(vel)), 'b-', 'linewidth', 3);  % show y = 0 on the grid
% set(gca,'fontsize',16);
% t1 = title(['f(x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v,t)  [DoGPack]']); 
% t1 = title(['f(t, x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v)']); 
% set(t1,'fontsize',16);

% axis([4.8 2*pi -1e-6 6e-6]);

% figure(1)

% % print a pretty picture!
% picname = strcat('Strong', outputdir);
% picname = strcat(picname, num2str(time) );
% print(1, '-depsc2', picname );

% picname = strcat(picname, 'slice');
% print(4, '-depsc2', picname );

figure(1);
