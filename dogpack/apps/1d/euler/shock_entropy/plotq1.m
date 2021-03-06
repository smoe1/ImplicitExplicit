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

% gas constant
%% ----
%% HACK
%% ----
gamma = 1.4;
%% ----

figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
%axis([-5 5 0 5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-5:2.5:5);
set(gca,'ytick',0:1:5);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(2);
clf;
press = (gamma-1).*(qsoln(:,5)-0.5*(qsoln(:,2).^2 + ...
                                    qsoln(:,3).^2 + qsoln(:, ...
                                                  4).^2)./qsoln(:,1));
pz=plot(xc,press,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
%axis([-5 5 0 12]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-5:2.5:5);
set(gca,'ytick',0:2:12);
set(gca,'fontsize',16);
t1 = title(['Pressure at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(xc,qsoln(:,2)./qsoln(:,1),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
%axis([-5 5 -0.5 3]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-5:2.5:5);
set(gca,'ytick',0:1:3);
set(gca,'fontsize',16);
t1 = title(['u^1(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
    
figure(4);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
%axis([1.0 3.2 2.5 5])
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-5:2.5:5);
set(gca,'ytick',0:0.5:5)
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
