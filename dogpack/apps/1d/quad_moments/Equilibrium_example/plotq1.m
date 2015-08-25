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

% sigma (width)
Zmx = length(qsoln(:,1));
alpha = zeros(Zmx,1);       
mu1 = zeros(Zmx,1);
mu2 = zeros(Zmx,1);
om1 = zeros(Zmx,1);
om2 = zeros(Zmx,1);

for mm=1:Zmx
    Zrho = rho(mm,1);
    Zp   = press(mm,1);
    Zu   = u(mm,1);
    Zq   = qheat(mm,1);
    Zr   = r(mm,1);
    Zalpha = GetAlpha(Zrho,Zp,Zq,Zr);
    alpha(mm,1) = Zalpha;
    if (alpha(mm,1)<1.0e-5)
        qt(mm,1) = 0.0;
    else
        qt(mm,1) = qheat(mm,1)/alpha(mm,1);
    end
    Zqt = qt(mm,1);
    mu1(mm,1) = u(mm,1) + 0.5*qt(mm,1)/press(mm,1) - sqrt(press(mm,1)*alpha(mm,1)/rho(mm,1) ...
                                              + 0.25*(qt(mm,1)/press(mm,1))^2);
    mu2(mm,1) = u(mm,1) + 0.5*qt(mm,1)/press(mm,1) + sqrt(press(mm,1)*alpha(mm,1)/rho(mm,1) ...
                                              + 0.25*(qt(mm,1)/press(mm,1))^2);;
    om1(mm,1) = rho(mm,1)*(u(mm,1)-mu2(mm,1))/(mu1(mm,1)-mu2(mm,1));
    om2(mm,1) = rho(mm,1)*(mu1(mm,1)-u(mm,1))/(mu1(mm,1)-mu2(mm,1));
    
    ZPC = [1, -2*Zqt/Zp, (8*Zrho*Zp^3*Zalpha+Zrho^2*Zqt^2-10*Zrho*Zp^3)/(Zrho^2*Zp^2),...
           (-10*Zrho*Zp^2*Zalpha*Zqt+12*Zrho*Zp^2*Zqt)/(Zrho^2*Zp^2),...
           (-3*Zrho*Zp*Zqt^2+3*Zrho*Zp*Zalpha*Zqt^2+15*Zp^4-24* ...
            Zp^4*Zalpha+10*Zp^4*Zalpha^2)/(Zrho^2*Zp^2), ...
           (-4*Zalpha^2*Zqt*Zp^3-6*Zqt*Zp^3+10*Zalpha*Zqt*Zp^3)/(Zrho^2*Zp^2)];
    ZF(mm,1:5) = sort(roots(ZPC));    
    
end

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

figure(6);
clf;
pz=plot(xc,alpha,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
%axis([0 1 0.47 0.63])
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',0.47:0.03:0.63);
set(gca,'fontsize',16);
t1 = title(['\alpha at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(7)
clf
pz=plot(xc,mu1,'b-');
set(pz,'linewidth',2);
set(pz,'markersize',4);
hold on;
pz=plot(xc,mu2,'r--');
set(pz,'linewidth',2);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
%axis([0 1 -1.5 1.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1.5:0.5:1.5);
set(gca,'fontsize',16);
t1 = title(['Abscissas at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(8)
clf
pz=plot(xc,om1,'b-');
set(pz,'linewidth',2);
set(pz,'markersize',4);
hold on;
pz=plot(xc,om2,'r--');
set(pz,'linewidth',2);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
%axis([0 1 0 1.25])
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.25:3);
set(gca,'fontsize',16);
t1 = title(['Weights at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(1)