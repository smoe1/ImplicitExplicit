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

%qex = exp(10*xc);
qex = 1+sin(2*pi*xc);

figure(1);
clf;
pt=plot(xc,qex,'r-');
set(pt,'linewidth',1.5);
hold on;
pz=plot(xc,qsoln,'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([0 1 -0.1 2.1]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',0:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

err = norm(qsoln-qex,2)/norm(qex,2);
disp(['   dx = ',num2str(dx,'%0.8e'),'         err = ',num2str(err,'%0.8e')]);
