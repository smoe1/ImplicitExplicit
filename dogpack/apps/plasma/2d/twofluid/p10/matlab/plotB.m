function plotB(mframe_in,outputdir_in)

  global mframe;
  if(nargin>=1)
    mframe=mframe_in;
  elseif(isempty(mframe))
    mframe=0;
  end

  global outputdir;
  if(nargin==2)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  
  format long e;
  
  %global mx my xlow xhigh ylow yhigh method;
  %method=zeros(1,7);
  %[meqn,mx,my,meth1,nplot,xl,yl,xc,yc]=readgrid(outputdir);
  %read_config_file([outputdir '/dogpack.data']);
  [nout,tfinal,dtv,cflv,nv,method,...
   meqn,mx,my,mbc,xlow,xhigh,ylow,yhigh,...
   mrestart,nstart,datafmt]...
   = read_dogpack_parameters2([outputdir '/dogpack.data']);
  [dx,dy,xl,yl,xc,yc]=get_grid(mx,my,xlow,xhigh,ylow,yhigh,method(1));

  basename = sprintf('%04d',mframe);
  fname = [outputdir, '/Q21:23.',  basename, '.dat'];
  [time,data]=read_from_file(fname,mx*method(1),my*method(1),3);
  clf;
  %hold off;
  hold on;
  %mx
  %my
  %method(1)
  %size(xc)
  %size(yc)
  %size(data)
  %error('stop here')
  plot_fieldlines(xc,yc,data(:,:,1),data(:,:,2),dx,dy);
  % represent the out-of-plane component with color
  plot_scalar(data(:,:,3),xl,yl,-1,-1,[-.03,.03]);
  % set skip so that roughly the same number of components are
  % displayed regardless of the size of the input
  plot_vector(xc,yc,data(:,:,1),data(:,:,2),-1,.4);
  hold off;
  set_axes(xlow,xhigh,ylow,yhigh);
  set_title('B',time);

  %B1 = data(:,:,1);
  %B2 = data(:,:,2);

  % calculate the accumulation integral of the flux
  % around the boundary to get the boundary
  % conditions for the Poisson problem
  % 
  % (will apply a linear shift to make the sum
  % come out to zero)
  %boundary_flx = B2(:,1);
  %
  %curlB_3=zeros
  %potential=fd2Poisson(
  %hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% garbage henceforth %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a stream function for the magnetic field on the interior
function phi = stream_func_from_vec_by_poisson(B1,B2,dx,dy)
  mx=size(B1,1);
  my=size(B1,2);
  nx=mx-2;
  ny=my-2;

  % need to check minus signs
  % - on each boundary of the cumulative boundary flux calculation
  % - in forcing function

  % create boundary data for poisson problem
  %
  bound_end=2*my+2*nx+1;
  diff_bound=zeros(bound_end,1);
  B1_left = B1(1,my:-1:1);
  B1_left = 0.5*(B1_left(1:end-1)+B1_left(2:end));
  diff_bound(2:my) = transpose(B1_left/dy);
  B2_bottom = B2(1:mx,1);
  B2_bottom = 0.5*(B2_bottom(1:end-1)+B2_bottom(2:end));
  diff_bound(my+1:mx+my-1) = B2_bottom/dx;
  B1_right = B1(mx,1:my);
  B1_right = 0.5*(B1_right(1:end-1)+B1_right(2:end));
  diff_bound(mx+my:nx+2*my) = B1_right/dy;
  B2_top = B2(mx:-1:1,my);
  B2_top = 0.5*(B2_top(1:end-1)+B2_top(2:end));
  diff_bound(nx+2*my+1:end) = B2_top/dx;
  boundary=cumsum(diff_bound);
  % hopefully the sum is (close to) zero, but in case not
  % we subtract a linear function so that it will be.
  boundary = boundary+((boundary(1)-boundary(end))/(bound_end-1))*transpose(0:bound_end-1);

  % define forcing function
  B2x = (B2(3:end,2:end-1)-B2(1:end-2,2:end-1))/(2*dx);
  B1y = (B1(2:end-1,3:end)-B1(2:end-1,1:end-2))/(2*dy);
  forcing = B2x-B1y;

  phi = fd2Poisson(forcing,dx,dy,boundary);

end
