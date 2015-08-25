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

% parameters
fids = fopen('param.data','r');
gamma      = fscanf(fids,'%g',1);  junk = fgets(fids);
mass_ratio = fscanf(fids,'%g',1);  junk = fgets(fids);
debye      = fscanf(fids,'%g',1);  junk = fgets(fids);
cs_light   = fscanf(fids,'%g',1);  junk = fgets(fids);
larmor_radius = fscanf(fids,'%g',1);  junk = fgets(fids);
fclose(fids);
	
% show B2 values
%
fig=2;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
pz=plot(xc,qsoln(:,17),'b-');
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',-2.0:0.4:2.0);
axis([0 1 -1.5 1.5]);
t1 = title(['B_2 at t = ',num2str(time),', r_g=', ...
           num2str(larmor_radius),'     [DoGPack]']);
    
% show B3 values
%
fig=3;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
pz=plot(xc,qsoln(:,18),'b-');
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',-2.0:0.4:2.0);
axis([0 1 -1.5 1.5]);
t1 = title(['B_3 at t = ',num2str(time),', r_g=',...
           num2str(larmor_radius),'     [DoGPack]']);


% compare densities
%
fig=4;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
%set(gca,'fontweight','bold');
set(gca,'linewidth',2);
pz=plot(xc,qsoln(:,1)+qsoln(:,11),'b-');
set(pz,'linewidth',2);
hold off;
axis on; box on; grid on;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',0:0.2:1.2);
axis([0 1 0 1.2]);
t1 = title(['Density at t = ',num2str(time),', r_g=',...
           num2str(larmor_radius),'     [DoGPack]']);

% compute pressures
%

P11_i = qsoln(:,5)  - qsoln(:,2).*qsoln(:,2)./qsoln(:,1);
P22_i = qsoln(:,8)  - qsoln(:,3).*qsoln(:,3)./qsoln(:,1);
P33_i = qsoln(:,10) - qsoln(:,4).*qsoln(:,4)./qsoln(:,1);
P12_i = qsoln(:,6)  - qsoln(:,2).*qsoln(:,3)./qsoln(:,1);
P13_i = qsoln(:,7)  - qsoln(:,2).*qsoln(:,4)./qsoln(:,1);
P23_i = qsoln(:,9)  - qsoln(:,3).*qsoln(:,4)./qsoln(:,1);

KE_e = (qsoln(:,12).*qsoln(:,12) ...
      + qsoln(:,13).*qsoln(:,13) ...
      + qsoln(:,14).*qsoln(:,14))./(qsoln(:,11)*2.);
press_e = (gamma-1)*(qsoln(:,15) - KE_e);

P11 = P11_i + press_e;
P12 = P12_i;
P13 = P13_i;
P22 = P22_i + press_e;
P23 = P23_i;
P33 = P33_i + press_e;

press_i = ( P11_i + P22_i + P33_i )/3.0;
press = press_i + press_e;
    
% compute eigenvalues of ion pressure tensor
%
a1=-P11_i-P22_i-P33_i;
a2=P11_i.*P22_i+P11_i.*P33_i+P22_i.*P33_i-P12_i.^2-P13_i.^2-P23_i.^2;
a3=P11_i.*P23_i.^2+P22_i.*P13_i.^2+P33_i.*P12_i.^2 ...
  -2*P12_i.*P13_i.*P23_i-P11_i.*P22_i.*P33_i;
Q=(3*a2-a1.^2)/9;
R=(9*a1.*a2-27*a3-2*a1.^3)/54;
discriminant=Q.^3+R.^2;
% discriminant should be negative
%find(discriminant>0)
root_disc = sqrt(discriminant); % should be imaginary
S=(R+root_disc).^(1/3);
% could just let T be the complex conjugate of S.
T=(R-root_disc).^(1/3);
% should be real
S_plus_T=real(S+T);
imag_S_minus_T = imag(S-T);
eig3i=S_plus_T-a1/3;
partA = -S_plus_T/2-a1/3;
partB = (sqrt(3)/2)*imag_S_minus_T;
eig2i= partA+partB;
eig1i= partA-partB;

% compute eigenvalues of total pressure tensor
%
a1=-P11-P22-P33;
a2=P11.*P22+P11.*P33+P22.*P33-P12.^2-P13.^2-P23.^2;
a3=P11.*P23.^2+P22.*P13.^2+P33.*P12.^2 ...
  -2*P12.*P13.*P23-P11.*P22.*P33;
Q=(3*a2-a1.^2)/9;
R=(9*a1.*a2-27*a3-2*a1.^3)/54;
discriminant=Q.^3+R.^2;
% discriminant should be negative
%find(discriminant>0)
root_disc = sqrt(discriminant); % should be imaginary
S=(R+root_disc).^(1/3);
% could just let T be the complex conjugate of S.
T=(R-root_disc).^(1/3);
% should be real
S_plus_T=real(S+T);
imag_S_minus_T = imag(S-T);
eig3=S_plus_T-a1/3;
partA = -S_plus_T/2-a1/3;
partB = (sqrt(3)/2)*imag_S_minus_T;
eig2= partA+partB;
eig1= partA-partB;

% compare total pressures
%
fig=5;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
pz=plot(xc,eig3,'b-');
set(pz,'linewidth',2);
pr=plot(xc,eig1,'r-');
set(pr,'linewidth',2);
pr=plot(xc,eig2,'k-');
set(pr,'linewidth',2);
pr=plot(xc,press,'g--');
set(pr,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',-0.0:0.4:2.0);
axis([0 1  0.0 1.5]);
t1 = title(['roots at t = ',num2str(time),'     [DoGPack]']);
t1 = title(['Total Pressure at t = ',...
           num2str(time),', r_g=',num2str(larmor_radius),'     [DoGPack]']);
legend('P''33','P''11','P''22','10-mom. p');

% compare ion pressures
%
fig=6;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
pz=plot(xc,eig3i,'b-');
set(pz,'linewidth',2);
pr=plot(xc,eig1i,'r-');
set(pr,'linewidth',2);
pr=plot(xc,eig2i,'k-');
set(pr,'linewidth',2);
pr=plot(xc,press_i,'g--');
set(pr,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',-0.0:0.2:2.0);
axis([0 1  0.0 1.0]);
t1 = title(['roots at t = ',num2str(time),'     [DoGPack]']);
t1 = title(['Ion Pressure at t = ',num2str(time),', r_g=',...
           num2str(larmor_radius),'     [DoGPack]']);
legend('P''33_i','P''11_i','P''22_i','10-mom. p_i');

% compare elc pressures
%
fig=7;
if(ishandle(fig))
  set(0,'CurrentFigure',fig);
else
  figure(fig);
end
clf;
hold on;
pz=plot(xc,press_e,'b-');
set(pz,'linewidth',2);
pr=plot(xc,press_e,'r-');
set(pr,'linewidth',2);
pr=plot(xc,press_e,'k-');
set(pr,'linewidth',2);
pr=plot(xc,press_e,'g--');
set(pr,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[2.0 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',-0.0:0.4:2.0);
axis([0 1  0.0 1.0]);
t1 = title(['roots at t = ',num2str(time),'     [DoGPack]']);
t1 = title(['Electron Pressure at t = ',num2str(time),', r_g=',...
           num2str(larmor_radius),'     [DoGPack]']);
legend('P''33_e','P''11_e','P''22_e','10-mom. p_e');

