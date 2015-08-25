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


figure(1);
clf;
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
pz=plot(xc,u,'bo');
set(pz,'markersize',6);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
%axis([0 1 -0.05 0.55]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['Velocity at t = ',num2str(time)]); 
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

zvec = zeros(length(qsoln(:,1)),5);
if (abs(time-0.15)<=1.0e-12)
   
    for mm=1:length(qsoln(:,1))
        Zrho = rho(mm,1);
        Zu   = u(mm,1);
        Zp   = press(mm,1);
        Zq   = qheat(mm,1);
        Zr   = r(mm,1);
        Zalpha = GetAlpha(Zrho,Zp,Zq,Zr );
        [z,NumIts,mfound] = NewtonSys(Zrho,Zu,Zp,Zq,Zalpha);
        
        zvec(mm,1:5) = z(1:5);        
        
        %error('stop in the name of love...');
    end
    
    figure(6);
    clf;
    p1=plot(xc,zvec(:,1),'b-');
    set(p1,'linewidth',1);
    hold on;
    p2=plot(xc,zvec(:,2),'r-');    
    set(p2,'linewidth',1);
    p3=plot(xc,zvec(:,3),'k--');    
    set(p3,'linewidth',1);
    p4=plot(xc,zvec(:,4),'r-');    
    set(p4,'linewidth',1);
    p5=plot(xc,zvec(:,5),'b-');    
    set(p5,'linewidth',1);
    hold off;
    
end

figure(1)