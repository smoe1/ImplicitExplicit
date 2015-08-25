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

figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,1),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off',...
           'contour','on','levels',0:0.002:0.04);
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-1.02 1.02 -1.02 1.02]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['\phi(x,y)']); 
set(t1,'fontsize',16);
caxis([0 0.04]);
hold off;


figure(2);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,2),'xystyle',shad,...
           'zdata',qsoln(:,2),'zstyle',cont,'colorbar','off','mesh','off',...
           'contour','on','levels',-0.25:0.025:0.25);
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-1.02 1.02 -1.02 1.02]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['E^1(x,y)']); 
set(t1,'fontsize',16);
caxis([-0.25 0.25]);
hold off;

figure(3);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,3),'xystyle',shad,...
           'zdata',qsoln(:,3),'zstyle',cont,'colorbar','off','mesh','off',...
           'contour','on','levels',-0.3:0.03:0.3);
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-1.02 1.02 -1.02 1.02]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['E^2(x,y)']); 
set(t1,'fontsize',16);
caxis([-0.3 0.3]);
hold off;

figure(1);