%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EXTRA FUNCTION, SPECIFIED BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

  % Add in the index for each node (vertex)
% for i=1:NumNodes
%   xn=node(i,1);
%   yn=node(i,2);
%     
%   t1=text(xn,yn,[num2str(i)]);
%   set(t1,'color',[0 0 1]);
%   set(t1,'fontweight','bold');
%   set(t1,'fontsize',14);
% end

 % Plot the centroid of each element
%for i=1:NumElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%end

 BndyList = [];    % List for saving each triangle that intersects curve of
                   % interest

 for i=1:NumElems

    % Compute the centroid and location of each vertex (node)
    x = node( tnode(i,:), 1 );   xmid = sum( x ) / 3;
    y = node( tnode(i,:), 2 );   ymid = sum( y ) / 3;

    % compute sgn of the distance from the vertices to the line.  E.g.,
    % determine if the poitns are 'above' or 'below' the line:
%   tmp    = curve_func( [x, y] ) ./ abs( curve_func( [x,y] ) );
    tmp    = sign( curve_func( [x,y] ) );
    I      = find( abs( curve_func( [x,y] ) ) < 1e-14 );
    tmp(I) = 0;

    % determine the 'distance' to the line.
    % triangles within '3' of the line contain part of the line on their
    % interior
    d = abs( sum( tmp ) );

    if( d < 3 )
     t1=text(xmid,ymid,[num2str(i)]);
     set(t1,'color',[0 0 1]);
     set(t1,'fontweight','bold');
     BndyList = [BndyList i];
    end

 end

% plot the intersected triangles:
% Can perform a contour fit on curve_func here if we really want ...
hold on;
xline = linspace( xmin, xmax );
plot( xline, xline.^2-0.5, '--k', 'Linewidth', 2 );
hold off;


%% -------------------------------------------------- %%
% search for a single point on the curve:
%% -------------------------------------------------- %%


xpt = 0.0;
pt  = [xpt, -0.5];

% plot the point we're searching for
hold on;
plot( [pt(1) pt(1)], [pt(2) pt(2)], 'ko', 'LineWidth', 3 );
hold off;

fprintf(1,'Searching for a single point (%f,%f) in the grid\n', pt(1),pt(2) );
NumBndTriangles = length( BndyList );
for i=1:NumBndTriangles

    % Compute the centroid and location of each vertex (node)
    x = node( tnode(BndyList(i),:), 1 );   xmid = sum( x ) / 3;
    y = node( tnode(BndyList(i),:), 2 );   ymid = sum( y ) / 3;

    dx = pt(1) - xmid;  dy = pt(2) - ymid;
    A  = [ [x(2)-x(1), x(3)-x(1)]; [y(2)-y(1), y(3)-y(1)] ];
    spts = A \ [dx,dy]';

    xi = spts(1);  eta = spts(2);

    % check if this point is inside the triangle:
    one_third = 1./3.;
    if( xi > -one_third & eta > -one_third & eta+xi < one_third )
        fprintf(1,'    found the pt (%f,%f) in triangle number %d\n', ...
        pt(1), pt(2), BndyList(i));
    end
end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;
%
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[0 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',16);
%end

%for i=1:NumPhysElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%end

%for j=1:NumGhostElems
%  i = NumPhysElems+j;
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  k = ghost_link(j);
%  
%  t2=text(xmid,ymid,[num2str(k)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;
%
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end

%figure(1)
%hold on;
%s=transpose(linspace(0,0.3,20));
%for i=1:NumEdges
%  
%  x1 = edge(i,1);
%  x2 = edge(i,3);
%  y1 = edge(i,2);
%  y2 = edge(i,4);
%  
%  xmid = (x1+x2)/2;
%  ymid = (y1+y2)/2;
%  
%  L = sqrt((x2-x1)^2 + (y2-y1)^2);
%  n1 = (y2-y1)/L;
%  n2 = (x1-x2)/L;
% 
%  xv=xmid+s*n1;
%  yv=ymid+s*n2;
%
%  pvq=plot(xv,yv,'k--');
%  set(pvq,'linewidth',2);
%  
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end
%hold off;
%


%for i=1:NumPhysElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%  set(t1,'fontsize',14);
%end

%for i=1:NumPhysNodes
%  t1=text(node(i,1),node(i,2),[num2str(i)]);
%  set(t1,'color',[1 0 1]);
%  set(t1,'fontweight','bold');
%  set(t1,'fontsize',14);
%end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;
%
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end

figure(1)
