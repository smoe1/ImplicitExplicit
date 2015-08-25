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

for mmm=1:meqn
    qlow(mmm)  = min(qsoln(:,mmm));
    qhigh(mmm) = max(qsoln(:,mmm));
    qeps(mmm)  = 0.015*(qhigh(mmm)-qlow(mmm));
end

figure(1);
clf;
p1=plot(xc,qsoln(:,1),'b-');
set(p1,'linewidth',2);
axis on; box on; grid off;
axis([xlow xhigh qlow(1)-qeps(1) qhigh(1)+qeps(1)]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['E(t,x) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(2);
clf;
p2=plot(xc, -qsoln(:,2),'r-');
set(p2,'linewidth',2);
axis on; box on; grid off;
axis([xlow xhigh -qhigh(2)-qeps(2) -qlow(2)+qeps(2)]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['-\phi(t,x) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
