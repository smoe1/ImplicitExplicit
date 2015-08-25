%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                         mx, my, mz:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh,zlow,zhigh]:  min/max values of grid
%%                               meqn:  number of equations
%%                               maux:  number of aux components
%%                              meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc,zc): grid points (cell centers), size = (mx,my,mz)
%%       (xl,yl,zl): grid points (lower left cell corners), size = (mx+1,my+1,mz+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,mz,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,mz,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,mz+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,mz+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
clf;
fv = isosurface(xc,yc,zc,qsoln,0.01);
p=patch(fv);
axis on; box off; grid on;
axis([0 1 0 1 0 1]);
isonormals(xc,yc,zc,qsoln,p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3);
%camlight;
%lighting gouraud;
l1=xlabel('x');
l2=ylabel('y');
l3=zlabel('z');
set(gca,'fontsize',16);
set(l1,'fontsize',16);
set(l2,'fontsize',16);
set(l3,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);
set(gca,'ztick',0:0.25:1);

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

qex = zeros(mx,my,mz);
for i=1:mx
    for j=1:my
        for k=1:mz
            xx = xc(i,j,k);
            yy = yc(i,j,k);
            zz = zc(i,j,k);
            qex(i,j,k) = qexfunc(time,xx,yy,zz,params);
        end
    end
end

err = norm(reshape(qex(1:mx,1:my,1:mz)-qsoln(1:mx,1:my,1:mz),mx*my*mz,1),2)/...
      norm(reshape(qex(1:mx,1:my,1:mz),mx*my*mz,1),2);

disp(' ');
disp([' Error:  ',num2str(err,'%10.5e')]);
disp(' ');

figure(1);