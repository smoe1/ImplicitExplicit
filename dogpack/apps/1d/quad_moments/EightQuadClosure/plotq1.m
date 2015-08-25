%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% desnity
rho = qsoln(:,1);

% velocity
u1 = qsoln(:,2)./qsoln(:,1);
u2 = qsoln(:,3)./qsoln(:,1);
u3 = qsoln(:,4)./qsoln(:,1);

% pressure 
P11 = qsoln(:,5)-(qsoln(:,2).^2)./qsoln(:,1);
P22 = qsoln(:,8)-(qsoln(:,3).^2)./qsoln(:,1);
P33 = qsoln(:,10)-(qsoln(:,4).^2)./qsoln(:,1);

% heat flux
Q111 = qsoln(:,11) - 3.0*P11.*u1 - rho.*u1.^3;
Q222 = qsoln(:,17) - 3.0*P22.*u2 - rho.*u2.^3;
Q333 = qsoln(:,20) - 3.0*P33.*u3 - rho.*u3.^3;



figure(1);
clf;
pz=plot(xc,rho,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 3.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.2:1);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(2);
clf;
pz=plot(xc,u1,'bo');
set(pz,'markersize',4);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([0 1 -1.0 1.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['Mean velocity u1 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(xc,u2,'bo');
set(pz,'markersize',4);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([0 1 -1.0 1.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['Mean velocity u2 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(4);
clf;
pz=plot(xc,u3,'bo');
set(pz,'markersize',4);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([0 1 -1.0 1.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['Mean velocity u3 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(5);
clf;
pz=plot(xc,P11,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 3.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Pressure P11 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(6);
clf;
pz=plot(xc,P22,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 3.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Pressure P22 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(7);
clf;
pz=plot(xc,P33,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 3.0]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Pressure P33 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(8);
clf;
pz=plot(xc,Q111,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 2.50]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Heat fl Q111 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(9);
clf;
pz=plot(xc,Q222,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 2.50]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Heat fl Q222 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(10);
clf;
pz=plot(xc,Q333,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.0 2.50]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Heat fl Q333 at t = ',num2str(time)]); 
set(t1,'fontsize',16);


%for nn=1:5
%figure(nn+9)
%clf
%pz=plot(xc,Zconvexity(:,nn),'r-');
%set(pz,'linewidth',2);
%set(gca,'plotboxaspectratio',[1.5 1 1]);
%end
