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


q1ex = -cos(sqrt(5)*pi*xc)*sqrt(5)*pi-(2/5)*sin((6/5)*sqrt(5)*pi)+...
           (2/5)*sin((37/10)*sqrt(5)*pi)+2;

xtmp = 3.7;
q1   = -cos(sqrt(5)*pi*xtmp)*sqrt(5)*pi-(2/5)*sin((6/5)*sqrt(5)*pi)+...
           (2/5)*sin((37/10)*sqrt(5)*pi)+2

q2   = (6/5)*cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi-sin((6/5)*sqrt(5)*pi)...
       -cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi*xtmp+sin(sqrt(5)*pi*xtmp)...
       +(cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi+(2/5)*sin((6/5)*...
       sqrt(5)*pi)-(2/5)*sin((37/10)*sqrt(5)*pi)-2)*(xtmp-6/5)+2

q2ex = (6/5)*cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi-sin((6/5)*sqrt(5)*pi)...
       -cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi*xc+sin(sqrt(5)*pi*xc)...
       +(cos((6/5)*sqrt(5)*pi)*sqrt(5)*pi+(2/5)*sin((6/5)*...
       sqrt(5)*pi)-(2/5)*sin((37/10)*sqrt(5)*pi)-2)*(xc-6/5)+2;

figure(1);
clf;
pt=plot(xc,q1ex,'r-');
set(pt,'linewidth',1.5);
hold on;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([1.2 3.7 -6 10]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',1:0.5:4);
set(gca,'ytick',-20:4:20);
set(gca,'fontsize',16);
t1 = title(['q^1(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(2);
clf;
pt=plot(xc,q2ex,'r-');
set(pt,'linewidth',1.5);
hold on;
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([1.2 3.7 -5 3]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',1:0.5:4);
set(gca,'ytick',-20:2:20);
set(gca,'fontsize',16);
t1 = title(['q^2(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(1);

dx = xc(2)-xc(1);
err1 = norm(q1ex-qsoln(:,1),2)/norm(q1ex,2);
err2 = norm(q2ex-qsoln(:,2),2)/norm(q2ex,2);

disp(' ');
disp(['   dx = ',num2str(dx,'%0.8e'),'         err1 = ', ...
      num2str(err1,'%0.8e')]);
disp(['   dx = ',num2str(dx,'%0.8e'),'         err2 = ', ...
      num2str(err2,'%0.8e')]);
disp(' ');
