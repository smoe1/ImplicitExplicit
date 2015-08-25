clear;

mx = 31;
my = 31;
mz = 31;

mxslice = 16;
myslice = 16;
mzslice = 16;

xlow  = 0;
xhigh = 1;
ylow  = 0;
yhigh = 1;
zlow  = 0;
zhigh = 1;

dx = (xhigh-xlow)/mx;
dy = (yhigh-ylow)/my;
dz = (zhigh-zlow)/mz;

xc = (xlow+dx/2):dx:(xhigh-dx/2);
yc = (ylow+dy/2):dy:(yhigh-dy/2);
zc = (zlow+dz/2):dz:(zhigh-dz/2);

[yc,xc,zc]=meshgrid(yc,xc,zc);

xeps = max([0.015*(xhigh-xlow),0.015*(yhigh-ylow),0.015*(xhigh-xlow),0.015*(zhigh-zlow)]);
yeps = xeps;
zeps = xeps;

q = zeros(mx,my,mz);
for i=1:mx
    for j=1:my
        for k=1:mz
            x = xc(i,j,k);
            y = yc(i,j,k);
            z = zc(i,j,k);
            r2 = (x-0.40)^2 + (y-0.50)^2 + (z-0.50)^2;
            r = sqrt(r2);
            if (r<0.3)
                q(i,j,k) = cos(5.0/3.0*pi*r)^6;
            end
        end
    end
end
  

figure(1);
clf;
fv = isosurface(xc,yc,zc,q,0.01);
p=patch(fv);
axis on; box off; grid on;
axis([0 1 0 1 0 1]);
isonormals(xc,yc,zc,q,p);
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

figure(2);
clf
xx1(1:mx,1:my) = xc(:,:,mzslice);
yy1(1:mx,1:my) = yc(:,:,mzslice);
qq1(1:mx,1:my) =  q(:,:,mzslice);
contourf(xx1,yy1,qq1,15,'k-');
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
hold on;
plot([xlow xhigh xhigh xlow xlow],[ylow ylow yhigh yhigh ylow],'k-');
hold off;
l1=xlabel('x');
l2=ylabel('y');
set(gca,'fontsize',16);
set(l1,'fontsize',16);
set(l2,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);

figure(3);
clf
xx2(1:mx,1:mz) = xc(:,myslice,:);
zz2(1:mx,1:mz) = zc(:,myslice,:);
qq2(1:mx,1:mz) =  q(:,myslice,:);
contourf(xx2,zz2,qq2,15,'k-');
axis('equal');
axis([xlow-xeps xhigh+xeps zlow-zeps zhigh+zeps]);
hold on;
plot([xlow xhigh xhigh xlow xlow],[zlow zlow zhigh zhigh zlow],'k-');
hold off;
l1=xlabel('x');
l2=ylabel('z');
set(gca,'fontsize',16);
set(l1,'fontsize',16);
set(l2,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);

figure(4);
clf
yy3(1:my,1:mz) = yc(mxslice,:,:);
zz3(1:my,1:mz) = zc(mxslice,:,:);
qq3(1:my,1:mz) =  q(mxslice,:,:);
contourf(yy3,zz3,qq3,15,'k-');
axis('equal');
axis([ylow-yeps yhigh+yeps zlow-zeps zhigh+zeps]);
hold on;
plot([ylow ylow yhigh yhigh ylow],[zlow zhigh zhigh zlow zlow],'k-');
hold off;
l1=xlabel('y');
l2=ylabel('z');
set(gca,'fontsize',16);
set(l1,'fontsize',16);
set(l2,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);

figure(1);
