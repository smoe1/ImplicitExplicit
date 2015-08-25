% plot field lines of divergenceless vector field
% by plotting contours of stream function
%
function plot_fieldlines(xc,yc,vecx,vecy,dx,dy)
  stream_func = stream_func_from_vec(vecx,vecy,dx,dy);
  %stream_func = stream_func-stream_func(1,1);
  plot_contour(xc,yc,stream_func);
end

function plot_contour(xc,yc,phi)
   cx = [-xc(end:-1:1,end:-1:1),-xc(end:-1:1,:);...
         xc(:,end:-1:1),xc];
   cy = [-yc(end:-1:1,end:-1:1), yc(end:-1:1,:);...
         -yc(:,end:-1:1),yc];
   full_phi = [phi(end:-1:1,end:-1:1),phi(end:-1:1,:);...
         phi(:,end:-1:1),phi];
   maxp = max(max(phi));
   minp = min(min(phi));
   contour_step=(maxp-minp)/31;
   upper=0:contour_step:maxp;
   lower=0:-contour_step:minp;
   levels=[lower(end:-1:2),upper];
   %levels=-15:15/15.;
   contour(cx,cy,full_phi,30,'-k');
   %contour(cx,cy,full_phi,levels,'-k');
   % should use a higher-order estimate of phi(0,0)
   contour(cx,cy,full_phi-phi(1,1),[0,1000],'--m');
end

function phi = stream_func_from_vec(B1,B2,dx,dy)
  mx=size(B1,1);
  my=size(B1,2);
  % define mesh representing integral of B2 from
  % each point's left neighbor to itself
  B2_cell_xint = dx*padarray(0.5*(B2(1:end-1,:)+B2(2:end,:)),[1,0],'pre');
  %B2_cell_xint = dx*0.5*(B2(1:end-1,:)+B2(2:end,:));
  % define mesh representing integral of B1 from
  % each point's upper neighbor to itself
  %B1_cell_yint = dy*padarray(0.5*(B1(:,1:end-1)+B1(:,2:end)),[0,1],'post');
  B1_cell_yint = dy*0.5*(B1(:,1:end-1)+B1(:,2:end));

  % for the sake of numerical stability we define the stream function to
  % be a weighted average of the path integral from the upper left to the given
  % point over all possible paths that move down and to the right along
  % lines between neighboring cell centers.  Inductively we say that the
  % value in each cell is an average of the value of the value in the
  % cell to the left plus the integral up to the present cell
  % with the value in the cell above plus the integral up to the present cell.
  %
  phi=zeros(size(B1));
  % first fill in the boundaries
  for ix=2:mx
    phi(ix,my) = phi(ix-1,my) + B2_cell_xint(ix,my);
  end
  for iy=my-1:-1:1
    phi(1,iy) = phi(1,iy+1) + B1_cell_yint(1,iy);
  end
  % now work down and to the right
  for iy=my-1:-1:1
  for ix=2:mx
    val_from_left = phi(ix-1,iy)+B2_cell_xint(ix,iy);
    val_from_abov = phi(ix,iy+1)+B1_cell_yint(ix,iy);
    phi(ix,iy) = 0.5*(val_from_left+val_from_abov);
  end
  end
end

