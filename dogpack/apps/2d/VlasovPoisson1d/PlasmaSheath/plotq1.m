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

for mmm=1:meqn
    qlow(mmm)  = min(qsoln(:,mmm));
    qhigh(mmm) = max(qsoln(:,mmm));
    qeps(mmm)  = 0.015*(qhigh(mmm)-qlow(mmm));
end

N = 1.0e13; % 1/m
L = 0.1; % m
T = 5.605424055e-9; % s
E0 = 1.809512620e4; % N/C
Phi0 = 1.809512620e3; % Nm/C

figure(1);
clf;
p1=plot(L*xc,0*ones(size(xc)),'k--');
set(p1,'linewidth',1);
hold on;
p2=plot(L*xc,E0*qsoln(:,1),'b-');
set(p2,'linewidth',2);
axis on; box on; grid off;
axis([0 0.1 -810 810]); 
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.02:0.1);
set(gca,'ytick',-800:200:800);
t1 = title(['E(t,x) at t = ',num2str(T*time)]); 
set(t1,'fontsize',16);
hold off;

h_mx = floor( mx/2 );
Phi_scale = Phi0*0.5*(qsoln(h_mx,2)+qsoln(h_mx+1,2));

figure(2);
clf;
p1=plot(L*xc,Phi_scale*ones(size(xc)),'k--');
set(p1,'linewidth',1);
hold on;
p2=plot(L*xc,Phi0*qsoln(:,2),'r-');
set(p2,'linewidth',2);
axis on; box on; grid off;
axis([0 0.1 -3.1 0.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.02:0.1);
set(gca,'ytick',-3:0.5:0);
t1 = title(['\phi(t,x) at t = ',num2str(T*time)]); 
set(t1,'fontsize',16);
hold off;

figure(3);
clf;
p1=plot(L*xc,N*ones(size(xc)),'k--');
set(p1,'linewidth',1);
hold on;
p2=plot(L*xc,N*qsoln(:,3),'m-');
set(p2,'linewidth',2);
axis on; box on; grid off;
axis([0.0 0.1 0.0 1.25e13]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.02:0.1);
set(gca,'ytick',0:0.25e13:1.25e13);
t1 = title(['n_e(t,x) at t = ',num2str(T*time)]); 
set(t1,'fontsize',16);
hold off;

disp([['Minimum value of number density = ', num2str( N*min( qsoln(:,3) ), '%2.5e' ) ]]);

  % print the pretty pictures!
  % frame number = n1
% fname = strcat( strcat('sheath-1d-efield', num2str(n1, '%03d') ), '.eps' );
% print(1,'-depsc', fname );

% fname = strcat( strcat('sheath-1d-density', num2str(n1, '%03d') ), '.eps' );
% print(3,'-depsc', fname );


