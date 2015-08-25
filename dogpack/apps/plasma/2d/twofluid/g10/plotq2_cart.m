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
surf( xl, yl, qaug(:,:,m));
hold on;
surf(-xl, yl, qaug(:,:,m));
surf( xl,-yl, qaug(:,:,m));
surf(-xl,-yl, qaug(:,:,m));
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-6:2:6);
set(gca,'ytick',-2:2:2);

figure(2);
clf;
contour( xc, yc, qsoln(:,:,m),15,'k');
hold on;
contour(-xc, yc, qsoln(:,:,m),15,'k');
contour( xc,-yc, qsoln(:,:,m),15,'k');
contour(-xc,-yc, qsoln(:,:,m),15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
set(gca,'xtick',-6:2:6);
set(gca,'ytick',-2:2:2);

figure(1);