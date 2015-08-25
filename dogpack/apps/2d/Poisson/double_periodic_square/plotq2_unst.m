%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%                NumElems:  number of elements
%%            NumPhysElems:  number of elements (excluding ghost elements)
%%           NumGhostElems:  number of ghost elements
%%                NumNodes:  number of nodes
%%            NumPhysNodes:  number of nodes (excluding exterior to domain)
%%             NumBndNodes:  number of nodes on boundary
%%                NumEdges:  number of edges
%% [xlow,xhigh,ylow,yhigh]:  bounding box
%%                    node:  list of nodes in mesh
%%                   tnode:  nodes attached to which element
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (NumPhysElem,meqn)
%%           aux:  aux components sampled on mesh, size = (NumPhysElem,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_ex = zeros(NumPhysElems*points_per_dir^2,3);
for i=1:(NumPhysElems*points_per_dir^2)
  xm = xmid(i);
  ym = ymid(i);

  q_ex(i,1) = sin(2*pi*xm)*sin(2*pi*ym);
  q_ex(i,2) = -2*pi*cos(2*pi*xm)*sin(2*pi*ym);
  q_ex(i,3) = -2*pi*cos(2*pi*ym)*sin(2*pi*xm);
end

figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,1),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','on','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['\phi(x,y)']); 
set(t1,'fontsize',16);
caxis([-1 1]);
hold off;


figure(2);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,2),'xystyle',shad,...
           'zdata',qsoln(:,2),'zstyle',cont,'colorbar','on','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['E^1(x,y)']); 
set(t1,'fontsize',16);
caxis([-2*pi 2*pi]);
hold off;

figure(3);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,3),'xystyle',shad,...
           'zdata',qsoln(:,3),'zstyle',cont,'colorbar','on','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['E^2(x,y)']); 
set(t1,'fontsize',16);
caxis([-2*pi 2*pi]);
hold off;

figure(1);

err1 = norm(q_ex(:,1) - qsoln(:,1),2)/norm(q_ex(:,1),2);
err2 = norm(q_ex(:,2) - qsoln(:,2),2)/norm(q_ex(:,2),2);
err3 = norm(q_ex(:,3) - qsoln(:,3),2)/norm(q_ex(:,3),2);


disp(' ');
disp(['   L2-error in phi: ',num2str(err1,'%10.5e')]);
disp(['   L2-error in E1:  ',num2str(err2,'%10.5e')]);
disp(['   L2-error in E2:  ',num2str(err3,'%10.5e')]);
disp(' ');