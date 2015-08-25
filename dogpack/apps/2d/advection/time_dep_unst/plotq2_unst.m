%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%                NumElems:  number of elements
%%            NumPhysElems:  number of elements (excluding ghost elements)
%%           NumGhostElems:  number of ghost elements
%%                NumNodes:  number of nodes
%%            NumPhysNodes:  number of nodes (excluding exterior to domain)
%%             NumBndNodes:  number of nodes on boundary
%%                NumEdges:  number of edges
%% [xlow,xhigh,ylow,yhigh]:  bounding box
%%                    node:  list of nodes in mesh
%%                   tnode:  nodes attached to which element
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (NumPhysElem,meqn)
%%           aux:  aux components sampled on mesh, size = (NumPhysElem,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -------------------------------------------------------------------------- %%
%% Scatter plot of the exact Solution
%% -------------------------------------------------------------------------- %%

% Center of the circle:
xstar = 0.0;
ystar = 0.0;
for i=1:(NumPhysElems*points_per_dir^2)
  xm = xmid(i);
  ym = ymid(i);

  rscat(i,1) = sqrt((xm-xstar)^2 + (ym-ystar)^2);
  qex(i,1) = (0.75+0.25*cos(2*pi*time)) * (0.5+0.5*cos(pi*rscat(i,1))); 
end

% Plots for the exact solution
[rscat,Itmp] = sort(rscat);
qex(:,1)     = qex(Itmp,1);


%% -------------------------------------------------------------------------- %%
%% Construct a slice of the solution at y = 0
%% -------------------------------------------------------------------------- %%


% Grab the slice index for this problem:
[BndyList, Psort, MidPt, LftPt, RghtPt] = FindSliceIndexUnst2( );
NumTriOnLine = length( BndyList ); 


% Sample the function on each triangle (this won't be a uniform sampling ... )
[q_data, T] = read_state2_unst(datafmt, outputdir, n1, 'q', ...
    NumElems, NumPhysElems, meqn, kmax, 1:meqn);      

local_pts_per_dir = points_per_dir;
rline = zeros( local_pts_per_dir*NumTriOnLine, 1 );
qline = zeros( local_pts_per_dir*NumTriOnLine, 1 );
index = 0;
for ii=1:NumTriOnLine

  % Find the coefficients for this triangle:
  TriNum = BndyList( Psort(ii) );
  qtmp   = zeros(1, meqn, kmax);
  qtmp(1, :, : ) = q_data( TriNum, : , : );

  % Center of current triangle:
  % REMEMBER TO USE OLD HERE: DivideUnst2Mesh modifies node and tnode!!!
  x_tmp = node_old( tnode_old(TriNum,:), 1 );   x_m = sum( x_tmp ) / 3;
  y_tmp = node_old( tnode_old(TriNum,:), 2 );   y_m = sum( y_tmp ) / 3;

  % Triangle boundary
  x1 = x_tmp(1);  y1 = y_tmp(1);
  x2 = x_tmp(2);  y2 = y_tmp(2);
  x3 = x_tmp(3);  y3 = y_tmp(3);
  A = [ (x2-x1), (x3-x1); (y2-y1), (y3-y1) ];

  % Canonical points we wish to sample:
  s2d_tmp = zeros( local_pts_per_dir, 2 );
  x = zeros( local_pts_per_dir, 1 );
  xl_tmp = LftPt( Psort(ii), 1 );  xr_tmp = RghtPt( Psort(ii), 1 );
  for mi=1:local_pts_per_dir
    x(mi)         = xl_tmp + (mi-0.5)*( xr_tmp-xl_tmp )/local_pts_per_dir;
    assert( x(mi) >= xl_tmp & x(mi) <= xr_tmp );
    s2d_tmp(mi,:) = A \ [ x(mi)-x_m; -y_m ];
  end

  % don't reuse LegVals here!  (It's used in the main cycle in plotdog2 ... )
  LegVals_tmp = GetUnst2Legendre(kmax, s2d_tmp);
  qvals       = sample_state2_unst(qtmp, meth1, kmax, LegVals_tmp);

  for mp=1:local_pts_per_dir
    index = index+1;
    rline(index) = x(mp);
    qline(index) = qvals(mp, 1 );
  end
end
clear q_data;

%% -------------------------------------------------------------------------- %%
%% Plots
%% -------------------------------------------------------------------------- %%

% Color plot of the entire solution
figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-1.01 1.01 -1.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
hold off;

% Errors (useful for convergence studies):
err = norm(qex-qsoln(Itmp,1),2)/norm(qex,2);
disp(' ');
disp(['  Error = ',num2str(err,'%0.5e')]);
disp(' ');

% Add in the single line, y=0 to the graph:
%   hold on;
%   xmin = -1;  xmax = 1.;
%   xline = linspace( xmin, xmax );
%   plot( xline, 0.*xline, '--k', 'Linewidth', 2 );
%   plot(  MidPt(Psort, 1), 0*MidPt(:,1), 'go', 'Linewidth', 3 );
%   hold off;

% Black and white contour plot
figure(2)
clf;
H2=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle','off',...
           'contour','on','levels',12,'colorbar','off');
for i=1:length(H2)
  set(H2(i),'linewidth',2*0.5);
  set(H2(i),'color','k');
end
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

% Scatter plot of the radial components
figure(3);
clf;
pz=plot(rscat,qsoln(Itmp,m),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold on;
pt=plot(rscat,qex(:,m),'r-');
set(pt,'linewidth',2);
hold off; 
axis on; box on; grid off;
axis([0 1.0 -0.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(t,r) at t = ',...
            num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

% Single plot of a slice through y = 0
figure(4);
clf;
pz = plot( rline, (0.75+0.25*cos(2*pi*time)) * (0.5+0.5*cos(pi*rline)), 'r-' ); 
set(pz, 'linewidth', 2);
axis on; box on; grid off;
axis([-1 1.0 -0.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(t,r) at t = ',...
                num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
hold on;
p2 = plot( rline, qline, 'bo' );

figure(1)
