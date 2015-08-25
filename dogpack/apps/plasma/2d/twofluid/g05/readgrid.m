function [meqn,mx,my,method_order,nplot,xl,yl,xc,yc,datafmt]=readgrid(outputdir);

  %if(1)
  [ nout    , ...
    tfinal  , ...
    dtv     , ...
    cflv    , ...
    nv      , ...
    method  , ...
    meqn    , ...
    mx      , ...
    my      , ...
    mbc     , ...
    xlow    , ...
    xhigh   , ...
    ylow    , ...
    yhigh   , ...
    mrestart, ...
    nstart  , ...
    datafmt ] ...
  = read_dogpack_parameters2([outputdir '/dogpack.data']);
  %
  nplot=nout;
  [dx,dy,xl,yl,xc,yc]=get_grid(mx,my,xlow,xhigh,ylow,yhigh,method(1));

% deprecating qhelp and grid files
%
% else
% % meqn and nplot
% fids = fopen([outputdir '/qhelp.dat'],'r');
% meqn = fscanf(fids,'%d',1);
% nplot = fscanf(fids,'%d',1);
% method_order = fscanf(fids,'%d',1);
% mx    = fscanf(fids,'%d',1);
% my    = fscanf(fids,'%d',1);
% datafmt = 1;
% datafmt = fscanf(fids,'%d',1);
% fclose(fids);
% mx_out = method_order*mx;
% my_out = method_order*my;
% 
% % grid
% 
% fids = fopen([outputdir '/gridx.dat'],'r');
% xl = fscanf(fids,'%e',[inf]);
% xl = (reshape(xl,mx_out+1,my_out+1));
% fclose(fids);
% 
% fids = fopen([outputdir '/gridy.dat'],'r');
% yl = fscanf(fids,'%e',[inf]);
% yl = (reshape(yl,mx_out+1,my_out+1));
% fclose(fids);
% 
% fids = fopen([outputdir '/gridxc.dat'],'r');
% xc_t = fscanf(fids,'%e',[inf]);
% xc = (reshape(xc_t,mx_out,my_out));
% fclose(fids);
% 
% fids = fopen([outputdir '/gridyc.dat'],'r');
% yc_t = fscanf(fids,'%e',[inf]);
% yc = (reshape(yc_t,mx_out,my_out));
% fclose(fids);
% end
end

