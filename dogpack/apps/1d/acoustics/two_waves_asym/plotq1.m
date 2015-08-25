%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if time==0
  qex = qsoln;
end

figure(1);
clf;
pt=plot(xc,qex(:,1),'r-');
set(pt,'linewidth',1.5);
hold on;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['q1 = q_t at t = ',num2str(time),'     [DoGPack]']);
axis([-1 1 -0.1 1.1]);
set(gca,'xtick', -1:.25:1);
set(gca,'ytick',0:.25:1);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

figure(2);
clf;  
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['q2 = q_x at t = ',num2str(time),'     [DoGPack]']);  
axis([-1 1 -1 1]);
set(gca,'xtick', -1:.25:1);
set(gca,'ytick',-1:0.25:1);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

%if ( abs(time-round(time) - 2.0) < 1.0e-10 && time > 0)
if ( abs(time- 1.0) < 1.0e-10 && time > 0)
  err1 = norm(qsoln(:,1)-qex(:,1),1)/norm(qex(:,1),1);
  err2 = norm(qsoln(:,1)-qex(:,1),2)/norm(qex(:,1),2);
  disp(' ');
%  disp(['   dx = ',num2str(dx,'%0.8e'),['         err1 = ' ...
%                          ''],num2str(err1,'%0.8e')]);
%  disp(' ');
  disp(['   dx = ',num2str(dx,'%0.8e'),['         err2 = ' ...
                      ''],num2str(err2,'%0.8e')]);
  disp(' ');

  %  plot the exact solution on top of computed solution
  figure(1);
  hold on;
  pz=plot(xc,qex(:,1),'-k');
  hold off;
end

%max(qsoln)
