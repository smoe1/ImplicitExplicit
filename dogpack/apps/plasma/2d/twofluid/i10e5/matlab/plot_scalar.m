
% change this to concatenate into one big array and plot in one
% command to avoid patchy-looking result from each quadrant
% having its own color map (in case color_range is [0,0]).
%
% xflip: -1 means reverse
% yflip: -1 means reverse
%
function plot_scalar(qvals_in,xl,yl,xflip,yflip,color_range);
    [mx,my]=size(qvals_in);
    if(color_range(2)>0)
      qvals=color_range(2)*tanh(qvals_in/color_range(2));
      %qvals=color_range(2)*sign(qvals_in).*sqrt(tanh(abs(qvals_in)/color_range(2)));
    else
      qvals=qvals_in;
    end

    % we expect size(xl)==size(yl)==[mx+1,my+1];
    qplot=zeros(size(qvals)*2+1);
    % should decide convention for not flipping (1?),
    % but this works as long as the user is consistent
    if(xflip==yflip)
      qplot(1:mx,1:my)=qvals(end:-1:1,end:-1:1);
    else
      qplot(1:mx,1:my)=-qvals(end:-1:1,end:-1:1);
    end
    if(xflip==-1)
      qplot(mx+1:end-1,1:my)=-qvals(:,end:-1:1);
    else
      qplot(mx+1:end-1,1:my)= qvals(:,end:-1:1);
    end
    if(yflip==-1)
      qplot(1:mx,my+1:end-1)=-qvals(end:-1:1,:);
    else
      qplot(1:mx,my+1:end-1)= qvals(end:-1:1,:);
    end
    qplot(mx+1:end-1,my+1:end-1)=qvals;

    % lower left coordinates
    lx=zeros(size(qvals)*2+1);
    ly=zeros(size(qvals)*2+1);
    lx(1:mx,1:my)=-xl(mx+1:-1:2,my+1:-1:2);
    ly(1:mx,1:my)=-yl(mx+1:-1:2,my+1:-1:2);
    lx(mx+1:end,1:my)= xl(:,my+1:-1:2);
    ly(mx+1:end,1:my)=-yl(:,my+1:-1:2);
    lx(1:mx,my+1:end)=-xl(mx+1:-1:2,:);
    ly(1:mx,my+1:end)= yl(mx+1:-1:2,:);
    lx(mx+1:end,my+1:end)= xl;
    ly(mx+1:end,my+1:end)= yl;

    pcolor(lx,ly,qplot);

%    % old code
%    qplot=zeros(size(qvals)+1);
%    qplot(1:end-1,1:end-1)=qvals;
%    pcolor(xl,yl,qplot);
%    shading flat;
%    hold on;
%    if(yflip==-1) pcolor( xl,-yl,-qplot);
%    else pcolor( xl,-yl,qplot);
%    end;
%    shading flat;
%    if(xflip==-1) pcolor(-xl, yl,-qplot);
%    else pcolor(-xl, yl,qplot);
%    end
%    shading flat;
%    pcolor(-xl,-yl,(xflip*yflip)*qplot);
     shading flat;
    colormap('jet');
    %low=min(qvals);
    %hgh=max(qvals);
    %mag=max([abs(low),hgh]);
    %mag=max(mag,.04);
    if(nargin>=6)
      caxis(color_range);
    end
    %else
      %caxis([0,mag]);
    %end
    c1=colorbar;
    set(c1,'fontsize',30);
    set(gca,'fontsize',30);
end

