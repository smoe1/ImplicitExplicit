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

% fourth moment
r = qsoln(:,5) - 4.0*qheat.*u - 6.0*press.*u.^2 - rho.*u.^4;

% sigma (width)
Zmx = length(qsoln(:,1));
alpha = zeros(Zmx,1);       
mu1 = zeros(Zmx,1);
mu2 = zeros(Zmx,1);
om1 = zeros(Zmx,1);
om2 = zeros(Zmx,1);
qt  = zeros(Zmx,1);
ZZ  = zeros(Zmx,5);
ZF  = zeros(Zmx,5);
ZRmat = zeros(Zmx,5,5);
ZLmat = zeros(Zmx,5,5);
Zchar = zeros(Zmx,5);
Zconvexity = zeros(Zmx,5);
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

    %[ZZtmp,NumIts,mfound] = NewtonSys(Zrho,u(mm,1),Zp,Zq,alpha(mm,1));
    %ZZ(mm,1:5) = sort(ZZtmp);
    
    %ZPC = [1,0,...
    %       (8*Zrho*Zp^3*Zalpha^3-10*Zrho*Zp^3*Zalpha^2)/(Zrho^2*Zp^2*Zalpha^2),0,...
    %       (10*Zp^4*Zalpha^4+15*Zp^4*Zalpha^2-24*Zp^4*Zalpha^3)/(Zrho^2*Zp^2*Zalpha^2),0];
    ZPC = [1, -2*Zqt/Zp, (8*Zrho*Zp^3*Zalpha+Zrho^2*Zqt^2-10*Zrho*Zp^3)/(Zrho^2*Zp^2),...
           (-10*Zrho*Zp^2*Zalpha*Zqt+12*Zrho*Zp^2*Zqt)/(Zrho^2*Zp^2),...
           (-3*Zrho*Zp*Zqt^2+3*Zrho*Zp*Zalpha*Zqt^2+15*Zp^4-24* ...
            Zp^4*Zalpha+10*Zp^4*Zalpha^2)/(Zrho^2*Zp^2), ...
           (-4*Zalpha^2*Zqt*Zp^3-6*Zqt*Zp^3+10*Zalpha*Zqt*Zp^3)/(Zrho^2*Zp^2)];
    ZF(mm,1:5) = sort(roots(ZPC));
    
    tmp = [1,1,1,1,1;
           ZF(mm,1),  ZF(mm,2),  ZF(mm,3),  ZF(mm,4),  ZF(mm,5);
           ZF(mm,1)^2,ZF(mm,2)^2,ZF(mm,3)^2,ZF(mm,4)^2,ZF(mm,5)^2;
           ZF(mm,1)^3,ZF(mm,2)^3,ZF(mm,3)^3,ZF(mm,4)^3,ZF(mm,5)^3;
           ZF(mm,1)^4,ZF(mm,2)^4,ZF(mm,3)^4,ZF(mm,4)^4,ZF(mm,5)^4];
    
    tmp_inv = inv(tmp);
    
    ZRmat(mm,:,:) = tmp;
    ZLmat(mm,:,:) = tmp_inv;
    
    char_tmp = tmp_inv*[qsoln(mm,1);qsoln(mm,2);qsoln(mm,3);qsoln(mm,4);qsoln(mm,5)];

    Zchar(mm,:) = char_tmp;
    
    DD = 0.0000005;
    ODD2 = 1/(2.0*DD);
    
    lp1 = sort(NewtonSys(Zrho+DD,Zu,Zp,Zq,Zr));
    lm1 = sort(NewtonSys(Zrho-DD,Zu,Zp,Zq,Zr));
    
    lp2 = sort(NewtonSys(Zrho,Zu+DD,Zp,Zq,Zr));
    lm2 = sort(NewtonSys(Zrho,Zu-DD,Zp,Zq,Zr));
    
    lp3 = sort(NewtonSys(Zrho,Zu,Zp+DD,Zq,Zr));
    lm3 = sort(NewtonSys(Zrho,Zu,Zp-DD,Zq,Zr));
    
    lp4 = sort(NewtonSys(Zrho,Zu,Zp,Zq+DD,Zr));
    lm4 = sort(NewtonSys(Zrho,Zu,Zp,Zq-DD,Zr));
    
    lp5 = sort(NewtonSys(Zrho,Zu,Zp,Zq,Zr+DD));
    lm5 = sort(NewtonSys(Zrho,Zu,Zp,Zq,Zr-DD));
    
    DPDU = [1, -Zu/Zrho, Zu^2, -Zu^3+(3*Zu*Zp)/Zrho, Zu^4+(4*Zu*Zq)/Zrho;
            0, 1/Zrho,-2*Zu, 3*Zu^2-(3*Zp)/Zrho, -4*Zu^3-(4*Zq)/Zrho;
            0, 0, 1, -3*Zu, 6*Zu^2;
            0, 0, 0, 1, -4*Zu;
            0, 0, 0, 0, 1];
    
    for nn=1:5
        
        DLDP = ODD2*[lp1(nn)-lm1(nn);
                     lp2(nn)-lm2(nn);
                     lp3(nn)-lm3(nn);
                     lp4(nn)-lm4(nn);
                     lp5(nn)-lm5(nn)];
        
        DLDU = DPDU*DLDP;
        
        Zconvexity(mm,nn) = tmp(1,nn)*DLDU(1) + tmp(2,nn)*DLDU(2) + ...
            tmp(3,nn)*DLDU(3) + tmp(4,nn)*DLDU(4) + tmp(5,nn)*DLDU(5);
    end
    
end

sigma = (press./rho).*(1-alpha);

figure(1);
clf;
pz=plot(xc,rho,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 0.9 1.9]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',1:0.2:2);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(2);
clf;
pz=plot(xc,u,'bo');
set(pz,'markersize',4);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
%axis([0 1 -0.05 0.55]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-1:0.25:1);
set(gca,'fontsize',16);
t1 = title(['Mean velocity at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(xc,press,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
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
set(pz,'markersize',4);
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
pz=plot(xc,qsoln(:,5),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
%axis([0 1 -2 0.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-3:0.5:3);
set(gca,'fontsize',16);
t1 = title(['M_4 at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(6);
clf;
pz=plot(xc,alpha,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',4);
hold off;
axis on; box on; grid off;
axis([0 1 -0.1 1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',0:0.2:1);
set(gca,'fontsize',16);
t1 = title(['\alpha at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(7)
clf
pz=plot(xc,mu1,'b-');
set(pz,'linewidth',2);
set(pz,'markersize',6);
hold on;
pz=plot(xc,mu2,'r--');
set(pz,'linewidth',2);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([0 1 -6 6]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-6:2:6);
set(gca,'fontsize',16);
t1 = title(['Abscissas at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(8)
clf
pz=plot(xc,om1,'b-');
set(pz,'linewidth',2);
set(pz,'markersize',6);
hold on;
pz=plot(xc,om2,'r--');
set(pz,'linewidth',2);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([0 1 -0.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-3:0.25:3);
set(gca,'fontsize',16);
t1 = title(['Weights at t = ',num2str(time)]); 
set(t1,'fontsize',16);

%figure(9)
%clf
%pz=plot(xc,ZF(:,1),'b-');
%set(pz,'linewidth',2);
%hold on;
%pz=plot(xc,ZF(:,2),'r--');
%set(pz,'linewidth',2);
%pz=plot(xc,ZF(:,3),'g-');
%set(pz,'linewidth',2);
%pz=plot(xc,ZF(:,4),'c--');
%set(pz,'linewidth',2);
%pz=plot(xc,ZF(:,5),'m-');
%set(pz,'linewidth',2);
%hold off;

for nn=1:5
figure(nn+9)
clf
pz=plot(xc,Zconvexity(:,nn),'r-');
set(pz,'linewidth',2);
set(gca,'plotboxaspectratio',[1.5 1 1]);
axis([0 1 -1 1]);
end

figure(1)