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

  %%% linux machines don't choke with contourf like they do with surf %%%
% surf(xl,yl,qaug(:,:,m));
  pcolor(xl,yl,qaug(:,:,m));

  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis([0 4*pi -2*pi 2*pi]);
  set(gca,'fontsize',16);
  t1 = title(['f(t,x,v) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  shading flat;
  colorbar;
  caxis auto;
% caxis([0,0.37]);
  set(gca,'fontsize',16);
  set(gca,'xtick',-6:2:6);
  set(gca,'ytick',-6:2:6);

  disp([['  minimum value of f is = ', num2str(min(min(qsoln)), '%2.8e')]])

% figure(2)
% clf;
% if mod(mx,2) == 0 % even number of cells printed!
%   xcenter_index = mx / 2;
% else
%   xcenter_index = (mx+1) / 2;
% end

% vel = yc(xcenter_index,:);
% qslice = qsoln(xcenter_index,:)';
% pp2=plot(vel, zeros(size(vel)), 'b-');  % show y = 0 on the grid
% set(pp2,'linewidth', 2.5);
% hold on;
% pp1=plot(vel, qslice,'r-');
% set(pp1,'linewidth', 2.5);
% axis([yl(1,1) yl(1,my) -0.05 0.45]);
% hold on;
% 
% set(gca,'fontsize',16);
% %t1 = title(['f(x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v,t) at t = ',num2str(time),'     [DoGPack]']); 
% t1 = title(['f(t, x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v) at t = ',num2str(time)]); 
% set(t1,'fontsize',16);


% %%%%%%%%%% Vertical Cross section that indicates soln goes negative %%%%%%%%
% xcut = xc(:,1);
% xcenter_index = max( find( xcut < 4.717 ) );
%
% figure(3);
% clf;
% vel = yc(xcenter_index,:);
% qslice = qsoln(xcenter_index,:)';
% pp2=plot(vel, zeros(size(vel)), 'b-');  % show y = 0 on the grid
% set(pp2,'linewidth', 2.5);
% hold on;
% pp1=plot(vel, qslice,'r-');
% set(pp1,'linewidth', 2.5);
% axis([yl(1,1) yl(1,my) -0.05 0.45]);
% hold on;
% 
% set(gca,'fontsize',16);
% t1 = title(['f(t, x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v) at t = ',num2str(time)]); 
% set(t1,'fontsize',16);

% figure(1);

  % Note that phi is stored in phi(:,2), and E = phi(:,1).
  phi = read_electric_field( n1, points_per_dir, point_type );

% parameters for the axes:
qlow  = min(min(phi) );
qhigh = max(max(phi) );
qeps  = 0.015*(qhigh-qlow);


  % Plot the electric potential
  figure(3);
  clf;
  p2 = plot( xc(:,1), -phi(:,2), 'r-' );
  set(p2,'linewidth',2);
  axis on; box on; grid off;
  axis([xlow xhigh -qhigh-qeps -qlow+qeps]);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'fontsize',16);
  t1 = title(['-\phi(t,x) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);

  % Plot the function: f_T( E ) := f( T, x, v ).
  % E.g., plot the graph { (E, f_T) : E = E(T,x,v), f_T = f(T,x,v) }
  %
  figure(4);
  clf;

  % compute energy at each point:
  Em = zeros( mx*my, 1 ); fp = zeros( mx*my, 1 );
  Ep = zeros( mx*my, 1 ); fm = zeros( mx*my, 1 );
  kp = 0; km = 0;
  for j=1:my
  for i=1:mx
    if( yc(i,j) > 0 )
      kp = kp+1;
      Ep( kp ) = 0.5*yc(i,j)^2 + phi(i,2);
      fp( kp ) = qsoln(i,j);
    else
      km = km+1;
      Em( km ) = 0.5*yc(i,j)^2 + phi(i,2);
      fm( km ) = qsoln(i,j);
    end
  end
  end

  plot( Ep(1:kp), fp(1:kp), 'go' );
  hold on;
  plot( Em(1:km), fm(1:km), 'r+' );

% axis([-0.01 0.015 0.197 0.204]);
  hold off;

  figure(1);

  %%%% Print the pretty pictures!
% picname1 = strcat(outputdir, ['q', num2str(n1,'%02d'), '.eps'] );
% picname2 = strcat(outputdir, '_slice', '.eps');
% picname3 = strcat(outputdir, '_slice2', '.eps');
% print(1, '-depsc2', picname1);
% print(2, '-depsc2', picname2);
% print(3, '-depsc2', picname3);
