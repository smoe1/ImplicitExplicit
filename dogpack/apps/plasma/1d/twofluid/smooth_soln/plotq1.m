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

% mass ratio
fidmr=fopen('mass_ratio.dat');
mass_ratio = fscanf(fidmr,'%e',1);
fclose(fidmr);

% Exact solution
nfunc = 1.0 + 0.5*sin(4.0*pi*(xc-time));
afunc = cos(4.0*pi*(xc-time));
bfunc = sin(2.0*pi*(xc-time));

rho_i_X   =  nfunc;
u1_i_X    =  1;
u2_i_X    =  0;
u3_i_X    =  0;
press_i_X =  1;
En_i_X    = press_i_X./(5/3-1) + 0.5*rho_i_X.*(u1_i_X.^2+u2_i_X.^2+u3_i_X.^2);

rho_e_X   =  nfunc/mass_ratio;
u1_e_X    =  1*ones(size(nfunc));
u2_e_X    =  0*ones(size(nfunc));
u3_e_X    =  0*ones(size(nfunc));
press_e_X =  0.2*ones(size(nfunc))/mass_ratio;
En_e_X    = press_e_X./(5/3-1) + 0.5*rho_e_X.*(u1_e_X.^2+u2_e_X.^2+u3_e_X.^2);

B1_X      =  1*ones(size(nfunc));
B2_X      =  afunc;
B3_X      =  bfunc;
E1_X      =  0*ones(size(nfunc));
E2_X      =  bfunc;
E3_X      = -afunc;	  

% Approximate solution
rho_i   = qsoln(:,1);
u1_i    = qsoln(:,2)./rho_i;
u2_i    = qsoln(:,3)./rho_i;
u3_i    = qsoln(:,4)./rho_i;
En_i    = qsoln(:,5);
press_i = (5/3-1).*(En_i - 0.5.*rho_i.*(u1_i.^2+u2_i.^2+u3_i.^2));

rho_e   = qsoln(:,6);
u1_e    = qsoln(:,7)./rho_e;
u2_e    = qsoln(:,8)./rho_e;
u3_e    = qsoln(:,9)./rho_e;
En_e    = qsoln(:,10);
press_e = (5/3-1).*(En_e - 0.5.*rho_e.*(u1_e.^2+u2_e.^2+u3_e.^2));

B1 = qsoln(:,11);
B2 = qsoln(:,12);
B3 = qsoln(:,13);

E1 = qsoln(:,14);
E2 = qsoln(:,15);
E3 = qsoln(:,16);

MZ=3;
NZ=2;

% ions
figure(1);
clf;
subplot(MZ,NZ,1);
plot(xc,rho_i,'bo');
hold on;
plot(xc,rho_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.5:0.25:1.5);
axis([0 1 0.4 1.6]);
t1 = title(['Ion density at t = ',num2str(time)]);

subplot(MZ,NZ,2);
plot(xc,press_i,'bo');
hold on;
plot(xc,press_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.8:0.1:1.2);
axis([0 1 0.8 1.2]);
t1 = title(['Ion pressure at t = ',num2str(time)]);

subplot(MZ,NZ,3);
plot(xc,u1_i,'bo');
hold on;
plot(xc,u1_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.8:0.1:1.2);
axis([0 1 0.8 1.2]);
t1 = title(['Ion x-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,4);
plot(xc,u2_i,'bo');
hold on;
plot(xc,u2_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-0.2:0.1:0.2);
axis([0 1 -0.2 0.2]);
t1 = title(['Ion y-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,5);
plot(xc,u3_i,'bo');
hold on;
plot(xc,u3_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-0.2:0.1:0.2);
axis([0 1 -0.2 0.2]);
t1 = title(['Ion z-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,6);
plot(xc,En_i,'bo');
hold on;
plot(xc,En_i_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',1.6:0.2:2.4);
axis([0 1 1.6 2.4]);
t1 = title(['Ion energy at t = ',num2str(time)]);

% electrons
figure(2);
clf;
subplot(MZ,NZ,1);
plot(xc,rho_e*mass_ratio,'bo');
hold on;
plot(xc,rho_e_X*mass_ratio,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.5:0.25:1.5);
axis([0 1 0.4 1.6]);
t1 = title(['(Electron density) X (mass ratio) at t = ',num2str(time)]);

subplot(MZ,NZ,2);
plot(xc,press_e*5*mass_ratio,'bo');
hold on;
plot(xc,press_e_X*5*mass_ratio,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.8:0.1:1.2);
axis([0 1 0.8 1.2]);
t1 = title(['(Electron pressure) X (5) X (mass ratio) at t = ',num2str(time)]);

subplot(MZ,NZ,3);
plot(xc,u1_e,'bo');
hold on;
plot(xc,u1_e_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.8:0.1:1.2);
axis([0 1 0.8 1.2]);
t1 = title(['Electron x-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,4);
plot(xc,u2_e,'bo');
hold on;
plot(xc,u2_e_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-0.2:0.1:0.2);
axis([0 1 -0.2 0.2]);
t1 = title(['Electron y-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,5);
plot(xc,u3_e,'bo');
hold on;
plot(xc,u3_e_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-0.2:0.1:0.2);
axis([0 1 -0.2 0.2]);
t1 = title(['Electron z-velocity at t = ',num2str(time)]);

subplot(MZ,NZ,6);
plot(xc,En_e,'bo');
hold on;
plot(xc,En_e_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',1e-4:1e-4:3e-4);
axis([0 1 1e-4 3e-4]);
t1 = title(['Electron energy at t = ',num2str(time)]);

% electromagnetics
figure(3);
clf;
subplot(MZ,NZ,1);
plot(xc,B1,'bo');
hold on;
plot(xc,B1_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0.8:0.1:1.2);
axis([0 1 0.8 1.2]);
t1 = title(['B1 at t = ',num2str(time)]);

subplot(MZ,NZ,2);
plot(xc,E1,'bo');
hold on;
plot(xc,E1_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-0.2:0.1:0.2);
axis([0 1 -0.2 0.2]);
t1 = title(['E1 at t = ',num2str(time)]);

subplot(MZ,NZ,3);
plot(xc,B2,'bo');
hold on;
plot(xc,B2_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.2:0.4:1.2);
axis([0 1 -1.2 1.2]);
t1 = title(['B2 at t = ',num2str(time)]);

subplot(MZ,NZ,4);
plot(xc,E2,'bo');
hold on;
plot(xc,E2_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.2:0.4:1.2);
axis([0 1 -1.2 1.2]);
t1 = title(['E2 at t = ',num2str(time)]);

subplot(MZ,NZ,5);
plot(xc,B3,'bo');
hold on;
plot(xc,B3_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.2:0.4:1.2);
axis([0 1 -1.2 1.2]);
t1 = title(['B3 at t = ',num2str(time)]);

subplot(MZ,NZ,6);
plot(xc,E3,'bo');
hold on;
plot(xc,E3_X,'r-');
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',12);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.2:0.4:1.2);
axis([0 1 -1.2 1.2]);
t1 = title(['E3 at t = ',num2str(time)]);