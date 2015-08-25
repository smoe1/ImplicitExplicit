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

% exact steady-state solution
fids = fopen('exact_soln/shllw.dat','r');
if fids~=-1
  tmp = fscanf(fids,'%e %e %e %e',[4,inf]);
  fclose(fids);
  x_ex = transpose(tmp(1,:));
  h_ex = transpose(tmp(2,:));
  b_ex = transpose(tmp(3,:));
  Fr_ex = transpose(tmp(4,:));
  clear tmp;
end

% bottom topography
bot = aux(:,1);

figure(1);
clf;
if (fids~=-1)
  pex=plot(x_ex,h_ex+b_ex,'r-');
  set(pex,'linewidth',1.5);
  hold on;
end
pz=plot(xc,qsoln(:,1)+bot,'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([0 1 0 1.2]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.2:1.2);
set(gca,'fontsize',16);
t1 = title(['h(x,t) + b(x)  at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
hold on;
fill(xc,bot,'g');
hold off;

if (time==10 && fids~=-1)  
  hi = interp1(x_ex,h_ex,xc,'spline');
  
  err = dx*norm(abs(hi-qsoln(:,1)),2);
end

figure(2);
clf;
pex=plot(xc,ones(size(qsoln(:,1))),'r-');
set(pex,'linewidth',1.5);
hold on;
u = qsoln(:,2)./qsoln(:,1);
Froude = u./sqrt(qsoln(:,1));
pz=plot(xc,Froude,'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['Froude number at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);


figure(3);
clf;
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['Momentum at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
