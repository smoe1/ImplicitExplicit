%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EXTRA FUNCTION, SPECIFIED BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold off;

figure(2) 
hold on;

%for i=1:NumElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;

%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 1 0]);
%  set(t2,'fontweight','bold');
%end

for i=1:NumPhysElems
  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;

  t1=text(xmid,ymid,[num2str(i)]);
  set(t1,'color',[0 0 1]);
  set(t1,'fontweight','bold');
end

for j=1:NumGhostElems
  i = NumPhysElems+j;
  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;

  k = ghost_link(j);
  
  t2=text(xmid,ymid,[num2str(k)]);
  set(t2,'color',[1 0 0]);
  set(t2,'fontweight','bold');
end


hold off;
figure(1)