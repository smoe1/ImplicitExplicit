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
u = qsoln(:,2)./qsoln(:,1);

% pressure
press = qsoln(:,3)-(qsoln(:,2).^2)./qsoln(:,1);

% heat flux
qheat = qsoln(:,4) - 3*press.*u - rho.*u.^3;
qstar = qheat./((press./rho).^(3/2));

% fourth moment
r = qsoln(:,5) - 4.0*qheat.*u - 6.0*press.*u.^2 - rho.*u.^4;


if (mx==100 && abs(time-0.2)<=1.0e-10)
   fileID = fopen('output_2000/ref2000.dat','r');
   A = fscanf(fileID, '%e %e %e', [3 inf]);
   fclose(fileID);
   A = A';
   
   xc_2000 = A(:,1);
   rho_2000 = A(:,2);
   M1_2000 = A(:,3);
end


figure(1);
clf;
if (mx==100 && abs(time-0.2)<=1.0e-10)
    pz=plot(xc_2000,rho_2000,'r-');
    set(pz,'linewidth',1);
    set(pz,'markersize',6);
    hold on;
end
pz=plot(xc,rho,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
%axis([0 1 0.8 3.2]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(2);
clf;
if (mx==100 && abs(time-0.2)<=1.0e-10)
    pz=plot(xc_2000,M1_2000,'r-');
    set(pz,'linewidth',1);
    set(pz,'markersize',6);
    hold on;
end
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',6);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
%axis([0 1 -0.05 0.55]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['M_1 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(xc,press,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
%axis([0 1 0.8 3.2]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Pressure at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(4);
clf;
pz=plot(xc,qheat,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
%axis([0 1 -2 0.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Heat flux at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(5);
clf;
pz=plot(xc,r,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
%axis([0 1 -2 0.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['Fourth moment at t = ',num2str(time)]); 
set(t1,'fontsize',16);

if (mx==2000 && abs(time-0.2)<=1.0e-10)
   x_out = xc;
   rho_out = qsoln(:,1);
   M1_out = qsoln(:,2);
   
   fileID = fopen('output_2000/ref2000.dat','w');
   for i=1:mx
       fprintf(fileID,'%16.8e %16.8e %16.8e\n',x_out(i),rho_out(i), ...
               M1_out(i));
   end
   fclose(fileID);
end



figure(1)