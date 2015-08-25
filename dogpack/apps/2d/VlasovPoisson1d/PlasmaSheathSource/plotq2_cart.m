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

N = 1.0e13; % 1/m
L = 0.1; % m
T = 5.605424055e-9; % s
V = L/T; % m/s
E0 = 1.809512620e4; % N/C
Phi0 = 1.809512620e3; % Nm/C
F = 5.605424055e5; % s/m^2
temp = 5.526350206e-4;

% electron mass and k is something:
me = 9.10938188e-31;
ke = 1.38064881e-23;

% maxwellian parameters:
thetat  = 1.160451897e4;
rhot    = 9.10938188e-18;
rho   = rhot / sqrt(2*pi*me*ke*thetat);
theta = 0.5*me/(ke*thetat);


  figure(1);
  clf;
  pcolor(L*xl,V*yl,F*qaug(:,:,m));
  colormap('jet');
  axis on; box on; grid off;
  axis([-0.0001*L 1.0001*L -0.2001*V 0.2001*V]);
  set(gca,'xtick',0:0.25*L:1*L);
  set(gca,'ytick',-0.2*V:0.1*V:0.2*V);  
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'fontsize',16);
  t1 = title(['f(t,x,v) at t = ',num2str(T*time)]); 
  set(t1,'fontsize',16);

  shading flat;
  caxis([0 17*F]);
  c1=colorbar;
  set(c1,'ytick',0:2.5*F:17*F);
  set(c1,'fontsize',16);

  disp([['  minimum value of f is = ', num2str(min(min(qsoln)), '%2.8e')]])

  figure(2)
  clf;
  if mod(mx,2) == 0 % even number of cells printed!
    xcenter_index = mx / 2;
  else
    xcenter_index = (mx+1) / 2;
  end

  vel = yc(xcenter_index,:);
  qslice = qsoln(xcenter_index,:)';
  
  p0=plot(vel, zeros(size(vel)), 'k--');  % show y = 0 on the grid
  set(p0,'linewidth',1);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  hold on;
  p1=plot(V*vel, F*qslice,'r-');
  set(p1,'linewidth',2);
  axis([-0.2*V 0.2*V -1*F 18*F]);
  set(gca,'xtick',-3e6:1e6:3e6);
  set(gca,'ytick',0:2e6:10e6);
  hold off;
  set(gca,'fontsize',16);
  t1 = title(['f(t, x = ',num2str(L*xc(xcenter_index,5), '%1.3f'), ',v) at t = ',...
              num2str(T*time)]); 
  set(t1,'fontsize',16);

  % subtract out a maxwellian:
  figure(3)
  clf;
  if mod(mx,2) == 0 % even number of cells printed!
    xcenter_index = mx / 2;
  else
    xcenter_index = (mx+1) / 2;
  end

  vel = yc(xcenter_index,:);
  qslice = qsoln(xcenter_index,:)';
  
  p0=plot(vel, zeros(size(vel)), 'k--');  % show y = 0 on the grid
  set(p0,'linewidth',1);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  hold on;

  %p1=plot(V*vel, (F*qslice) - (exp(-0.5*(vel').^2 / T ) / sqrt(2*pi*T)),'r-');
  p1=plot(V*vel, (F*qslice) - (rho*exp( -theta * (V*vel').^2 )), 'r-');
  hold on;
  set(p1,'linewidth',2);
  axis([-0.2*V 0.2*V -1*F 1.0*F]);
  set(gca,'xtick',-3e6:1e6:3e6);
  set(gca,'ytick',-F:(2*F/6):F);
  plot( V*vel, zeros(size(vel)), '--k' );
  hold off;
  set(gca,'fontsize',16);
  t1 = title(['f(t, x = ',num2str(L*xc(xcenter_index,5), '%1.3f'), ',v) - f(t=0,v) at t = ',...
              num2str(T*time)]); 
  set(t1,'fontsize',16);

  % print the pretty pictures!
  % frame number = n1
  fname = strcat( strcat('sheath-1d-phase', num2str(n1, '%03d') ), '.eps' );
  print(1,'-depsc', fname );

  fname = strcat( strcat('sheath-1d-slice', num2str(n1, '%03d') ), '.eps' );
  print(3,'-depsc', fname );

