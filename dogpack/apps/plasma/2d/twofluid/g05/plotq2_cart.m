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
surf( xl, yl, qaug(:,:,1));
hold on;
surf(-xl, yl, qaug(:,:,1));
surf( xl,-yl, qaug(:,:,1));
surf(-xl,-yl, qaug(:,:,1));
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Ion density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

figure(2);
clf;
contour( xc, yc, qsoln(:,:,1),15,'k');
hold on;
contour(-xc, yc, qsoln(:,:,1),15,'k');
contour( xc,-yc, qsoln(:,:,1),15,'k');
contour(-xc,-yc, qsoln(:,:,1),15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Ion density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

figure(3);
clf;
surf( xl, yl, qaug(:,:,6));
hold on;
surf(-xl, yl, qaug(:,:,6));
surf( xl,-yl, qaug(:,:,6));
surf(-xl,-yl, qaug(:,:,6));
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Electron density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

figure(4);
clf;
contour( xc, yc, qsoln(:,:,6),15,'k');
hold on;
contour(-xc, yc, qsoln(:,:,6),15,'k');
contour( xc,-yc, qsoln(:,:,6),15,'k');
contour(-xc,-yc, qsoln(:,:,6),15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Electron density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

Viaug = sqrt( (qaug(:,:,2)./qaug(:,:,1)).^2 + (qaug(:,:,3)./qaug(:,:,1)).^2 + (qaug(:,:,4)./qaug(:,:,1)).^2 );
figure(5);
clf;
surf( xl, yl, Viaug);
hold on;
surf(-xl, yl, Viaug);
surf( xl,-yl, Viaug);
surf(-xl,-yl, Viaug);
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Ion |V_i| at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

Vi = sqrt( (qsoln(:,:,2)./qsoln(:,:,1)).^2 + (qsoln(:,:,3)./qsoln(:,:,1)).^2 + (qsoln(:,:,4)./qsoln(:,:,1)).^2 );
figure(6);
clf;
contour( xc, yc, Vi,15,'k');
hold on;
contour(-xc, yc, Vi,15,'k');
contour( xc,-yc, Vi,15,'k');
contour(-xc,-yc, Vi,15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Ion |V_i| at  at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

Veaug = sqrt( (qaug(:,:,7)./qaug(:,:,6)).^2 + (qaug(:,:,8)./qaug(:,:,6)).^2 + (qaug(:,:,9)./qaug(:,:,6)).^2 );
figure(7);
clf;
surf( xl, yl, Veaug);
hold on;
surf(-xl, yl, Veaug);
surf( xl,-yl, Veaug);
surf(-xl,-yl, Veaug);
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Electron |V_e| at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

Ve = sqrt( (qsoln(:,:,7)./qsoln(:,:,6)).^2 + (qsoln(:,:,8)./qsoln(:,:,6)).^2 + (qsoln(:,:,9)./qsoln(:,:,6)).^2 );
figure(8);
clf;
contour( xc, yc, Ve,15,'k');
hold on;
contour(-xc, yc, Ve,15,'k');
contour( xc,-yc, Ve,15,'k');
contour(-xc,-yc, Ve,15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Electron |V_e| at  at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

Baug = sqrt( qaug(:,:,11).^2 + qaug(:,:,12).^2 + qaug(:,:,13).^2 );
figure(9);
clf;
surf( xl, yl, Baug);
hold on;
surf(-xl, yl, Baug);
surf( xl,-yl, Baug);
surf(-xl,-yl, Baug);
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['|B| at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

B = sqrt( qsoln(:,:,11).^2 + qsoln(:,:,12).^2 + qsoln(:,:,13).^2 );
figure(10);
clf;
contour( xc, yc, B,15,'k');
hold on;
contour(-xc, yc, B,15,'k');
contour( xc,-yc, B,15,'k');
contour(-xc,-yc, B,15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['|B| at  at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);


Eaug = sqrt( qaug(:,:,14).^2 + qaug(:,:,15).^2 + qaug(:,:,16).^2 );
figure(11);
clf;
surf( xl, yl, Eaug);
hold on;
surf(-xl, yl, Eaug);
surf( xl,-yl, Eaug);
surf(-xl,-yl, Eaug);
yrbcolormap;
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['|E| at t = ',num2str(time)]); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);

E = sqrt( qsoln(:,:,14).^2 + qsoln(:,:,15).^2 + qsoln(:,:,16).^2 );
figure(12);
clf;
contour( xc, yc, E,15,'k');
hold on;
contour(-xc, yc, E,15,'k');
contour( xc,-yc, E,15,'k');
contour(-xc,-yc, E,15,'k');
plot([-xhigh xhigh xhigh -xhigh -xhigh],[-yhigh -yhigh yhigh yhigh -yhigh],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([-xhigh-xeps xhigh+xeps -yhigh-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['|E| at  at t = ',num2str(time)]); 
set(t1,'fontsize',16);
set(gca,'xtick',-12:3:12);
set(gca,'ytick',-6:3:6);


figure(1);

%mxhalf = floor(mx/2)+1;
%myhalf = floor(my/2)+1;

%figure(3)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,1),'b-');
%title('rho_i x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,1),'b-');
%title('rho_i y-slice');

%figure(4)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,4)./qsoln(:,myhalf,1),'b-');
%title('u3_i x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,4)./qsoln(mxhalf,:,1),'b-');
%title('u3_i y-slice');

%figure(5)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,6),'b-');
%title('rho_e x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,6),'b-');
%title('rho_e y-slice');

%figure(6)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,9)./qsoln(:,myhalf,6),'b-');
%title('u3_e x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,9)./qsoln(mxhalf,:,6),'b-');
%title('u3_e y-slice');

%figure(7)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,11),'b-');
%title('B1 x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,11),'b-');
%title('B1 y-slice');

%figure(8)
%subplot(2,1,1);
%plot(xc,qsoln(:,myhalf,12),'b-');
%title('B2 x-slice');

%subplot(2,1,2);
%plot(yc,qsoln(mxhalf,:,12),'b-');
%title('B2 y-slice');
