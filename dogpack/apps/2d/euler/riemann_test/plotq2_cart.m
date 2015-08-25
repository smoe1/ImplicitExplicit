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

% gas constant
fids  = fopen([outputdir,'/eulerhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/eulerhelp.dat  not found.']);
end
gamma_gas  = fscanf(fids,'%e',1);
OPT        = fscanf(fids,'%d',1);
fclose(fids);

figure(1);
clf;
pcolor(xl,yl,qaug(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

% This routine was stripped from plotdog1.m.  If the format of qhelp.dat
% changes, then this will need to be modified as well.
fids  = fopen(['../../../1d/euler/shock_tube/output/qhelp.dat'],'r');
if (fids>0)

    ndims1d    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    GridType1d = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);

    % more information from qhelp.dat.  This is particular to the 1D solver.
    meqn_1d    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    maux_1d    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    nplot_1d   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    kmax1d     = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    datafmt_1d = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    %
    mx_1d      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    xlow_1d    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    xhigh_1d   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    fclose(fids);

    mx_old_1d = mx_1d;
    mx_1d     = mx_1d*points_per_dir;
    dx_1d     = (xhigh_1d-xlow_1d)/mx_1d;

    % Linearly spaced points
    xc_1d     = transpose(linspace(xlow_1d + 0.5*dx_1d, xhigh_1d - 0.5*dx_1d, mx_1d ));
    s1d       = -1.0 + (2.0*(1:points_per_dir)-1.0)/points_per_dir;

    % phi( 1:points_per_dir, 1:kmax1d )
    phi_1d = GetCart1Legendre(meth1, s1d);

    fname = ['../../../1d/euler/shock_tube/output/',...
           num2str(n1+10000),'.dat'];
    fname(37) = 'q';
    fids = fopen(fname,'r');
    time_1d = fscanf(fids,'%e',1);
    qtmp = fscanf(fids,'%e',[1,inf]);
    fclose(fids);

    qtmp      = transpose(qtmp);
    qtmp      = reshape(qtmp, mx_old_1d, meqn_1d, kmax1d);
    qsoln_1d  = zeros(mx_1d, meqn_1d);

    for i=1:mx_old_1d
    for me=1:meqn_1d
    for ii=1:points_per_dir
        v1(1:kmax1d,1) = phi_1d(ii,:);
        v2(1:kmax1d,1) = qtmp(i,me,:);
        qsoln_1d((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
    end
    end
    end
    clear qtmp;
  
    figure(2);
    clf;
    if (OPT==1)
        pz=plot(reshape(xc,mx*my,1),reshape(qsoln(:,:,1),mx*my,1),'bo');
    else
        pz=plot(reshape(yc,mx*my,1),reshape(qsoln(:,:,1),mx*my,1),'bo');
    end
    set(pz,'linewidth',2);
    set(pz,'markersize',8);
    hold on;
    pr = plot(xc_1d,qsoln_1d(:,1),'r-');
    set(pr,'linewidth',2);
    hold off;
    axis on; box on; grid off;
    set(gca,'plotboxaspectratio',[1.5 1 1]);
    set(gca,'fontsize',16);
    set(gca,'xtick',0:0.25:1);
    set(gca,'ytick',0:0.5:3.5);
    axis([0 1 0 3.5]);
    t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']);
end

figure(1)
