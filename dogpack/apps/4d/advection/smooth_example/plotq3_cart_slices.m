%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%               mx, my, mz:  # of grid points
%%               dx, dy, dz:  grid spacing
%%   mx_old, my_old, mz_old:  # of grid points (original mesh before sub-dividing)
%%   dx_old, dy_old, dz_old:  grid spacing (original mesh before sub-dividing)
%%                     meqn:  number of equations
%%                     maux:  number of aux components
%%                    meth1:  spatial order of accuracy
%%                     time:  time of current slice
%%                        m:  solution component that user asked to plot
%% [xlow,xhigh,ylow,yhigh,zlow,zhigh]:  min/max values of grid
%%
%%   Slice information:
%%        numslices:  number of slices
%%        slicekind:  type of slice (1:z=const, 2:y=const, 3:x=const),
%%                    size = (numslices)
%%         sliceidx:  index of slice, size = (numslices)
%%              ms1:  max value of 1st index in slice, size = (numslices)
%%              ms2:  max value of 2nd index in slice, size = (numslices)
%%
%%   Grid information:
%%       (xc,yc,zc): grid points (cell centers), size = (mx,my,mz)
%%       (xl,yl,zl): grid points (lower left cell corners), size = (mx+1,my+1,mz+1)
%%       (xc_z, yc_z): grid points (cell centers) at constant z, size = (mx,my)
%%       (xc_y, zc_y): grid points (cell centers) at constant y, size = (mx,mz)
%%       (yc_x, zc_x): grid points (cell centers) at constant x, size = (my,mz)
%%       (xl_z, yl_z): grid points (lower left corners) at constant z, size = (mx+1,my+1)
%%       (xl_y, zl_y): grid points (lower left corners) at constant y, size = (mx+1,mz+1)
%%       (yl_x, zl_x): grid points (lower left corners) at constant x, size = (my+1,mz+1)
%%
%%   Solution information:
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mmax+1,mmax+1,meqn,numslices),
%%                 where mmax = max([mx,my,mz])
%%
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mmax+1,mmax+1,meqn,numslices),
%%                 where mmax = max([mx,my,mz])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ns=1:numslices

    figure(ns);
    clf;
    switch(slicekind(ns))
      case 1
        zconst = zlow + (sliceidx(ns)-0.5)*dz_old;
        qtmp1(1:ms1(ns)+1,1:ms2(ns)+1)  = qaug(1:ms1(ns)+1,1:ms2(ns)+1,m,ns);
        p1=pcolor(xl_z,yl_z,qtmp1);
        set(gca,'fontsize',16);
        t1=title(['Slice of q(',num2str(m),') at z = ',num2str(zconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('x');
        ylabel('y');
        axis('equal');
        axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
        shading flat;
        hold on;
        p2=plot([xlow xhigh xhigh xlow xlow],...
                [ylow ylow yhigh yhigh ylow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);
      case 2
        yconst = ylow + (sliceidx(ns)-0.5)*dy_old;
        qtmp2(1:ms1(ns)+1,1:ms2(ns)+1)  = qaug(1:ms1(ns)+1,1:ms2(ns)+1,m,ns);
        p1=pcolor(xl_y,zl_y,qtmp2);
        set(gca,'fontsize',16);
        t1=title(['Slice of q(',num2str(m),') at y = ',num2str(yconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('x');
        ylabel('z');
        axis('equal');
        axis([xlow-xeps xhigh+xeps zlow-zeps zhigh+zeps]);
        shading flat;
        hold on;
        p2=plot([xlow xhigh xhigh xlow xlow],...
                [zlow zlow zhigh zhigh zlow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);
      case 3
        xconst = xlow + (sliceidx(ns)-0.5)*dx_old;
        qtmp3(1:ms1(ns)+1,1:ms2(ns)+1)  = qaug(1:ms1(ns)+1,1:ms2(ns)+1,m,ns);
        p1=pcolor(yl_x,zl_x,qtmp3);
        set(gca,'fontsize',16);
        t1=title(['Slice of q(',num2str(m),') at x = ',num2str(xconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('y');
        ylabel('z');
        axis('equal');
        axis([ylow-yeps yhigh+yeps zlow-zeps zhigh+zeps]);
        shading flat;
        hold on;
        p2=plot([ylow yhigh yhigh ylow ylow],...
                [zlow zlow zhigh zhigh zlow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);
    end
    
end

fids = fopen([outputdir,'/qhelp_InitialParams.dat'],'r');
x0    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
y0    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
z0    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
width = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
u     = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
v     = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
w     = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
fclose(fids);
half_pi_owidth = 0.5*pi/width;
params = [x0,y0,z0,width,half_pi_owidth,u,v,w,...
          xlow,xhigh,ylow,yhigh,zlow,zhigh];

for ns=1:numslices

    figure(ns+numslices);
    clf;
    switch(slicekind(ns))
      case 1
        zconst = zlow + (sliceidx(ns)-0.5)*dz_old;      
        qex1 = zeros(ms1(ns)+1,ms2(ns)+1);
        qex1(1:ms1(ns),1:ms2(ns)) = qexact_slices(time,...
                                                  xc_z,...
                                                  yc_z,...
                                                  zconst,...
                                                  1,...
                                                  params);
        p1=pcolor(xl_z,yl_z,qex1);
        set(gca,'fontsize',16);
        t1=title(['Slice of exact solution at z = ',num2str(zconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('x');
        ylabel('y');
        axis('equal');
        axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
        shading flat;
        hold on;
        p2=plot([xlow xhigh xhigh xlow xlow],...
                [ylow ylow yhigh yhigh ylow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);
      case 2
        yconst = ylow + (sliceidx(ns)-0.5)*dy_old;
        qex2 = zeros(ms1(ns)+1,ms2(ns)+1);
        qex2(1:ms1(ns),1:ms2(ns)) = qexact_slices(time,...
                                                  xc_y,...
                                                  yconst,...
                                                  zc_y,...
                                                  2,...
                                                  params);
        p1=pcolor(xl_y,zl_y,qex2);
        set(gca,'fontsize',16);
        t1=title(['Slice of exact solution at y = ',num2str(yconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('x');
        ylabel('z');
        axis('equal');
        axis([xlow-xeps xhigh+xeps zlow-zeps zhigh+zeps]);
        shading flat;
        hold on;
        p2=plot([xlow xhigh xhigh xlow xlow],...
                [zlow zlow zhigh zhigh zlow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);        
      case 3
        xconst = xlow + (sliceidx(ns)-0.5)*dx_old;
        qex3 = zeros(ms1(ns)+1,ms2(ns)+1);
        qex3(1:ms1(ns),1:ms2(ns)) = qexact_slices(time,...
                                                  xconst,...
                                                  yc_x,...
                                                  zc_x,...
                                                  3,...
                                                  params);
        p1=pcolor(yl_x,zl_x,qex3);
        set(gca,'fontsize',16);
        t1=title(['Slice of exact solution at x = ',num2str(xconst),...
                  ' at time t = ',num2str(time)]);
        xlabel('y');
        ylabel('z');
        axis('equal');
        axis([ylow-yeps yhigh+yeps zlow-zeps zhigh+zeps]);
        shading flat;
        hold on;
        p2=plot([ylow yhigh yhigh ylow ylow],...
                [zlow zlow zhigh zhigh zlow],'k-');
        set(p2,'linewidth',2);
        hold off;
        c1=colorbar;
        set(c1,'fontsize',16);
    end
    
end

%% errors
err = [0;0;0];

kk=1;
err(kk) = norm(reshape(qex1(1:ms1(kk),1:ms2(kk))-qtmp1(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2)/...
          norm(reshape(qex1(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2);

kk=2;
err(kk) = norm(reshape(qex2(1:ms1(kk),1:ms2(kk))-qtmp2(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2)/...
          norm(reshape(qex2(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2);

kk=3;
err(kk) = norm(reshape(qex3(1:ms1(kk),1:ms2(kk))-qtmp3(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2)/...
          norm(reshape(qex3(1:ms1(kk),1:ms2(kk)),ms1(kk)*ms2(kk),1),2);

disp(' ');
kk=1;
disp([' Error in slice 1:  ',num2str(err(1),'%10.5e')]);
kk=2;
disp([' Error in slice 2:  ',num2str(err(2),'%10.5e')]);
kk=3;
disp([' Error in slice 3:  ',num2str(err(3),'%10.5e')]);
disp(' ');

figure(1);