% generic routine to plot data
%
function plotq(outputdir_in)
  global outputdir
  if(nargin>=1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  
  format long e;
  
  % TITLE INFO --------------------------------------------------------------
  disp(' ')
  disp('  DOGPACK 2D GRID PLOTTING ROUTINE')
  disp(' ')
  % END TITLE INFO ----------------------------------------------------------
  
  %%% read in data
  
  [mx,my,xlow,xhigh,ylow,yhigh,method1,meqn,nout,mx_out,my_out,datafmt]=get_plot_params2(outputdir);
  [dx,dy,xl,yl]=get_left_grid2(mx,my,xlow,xhigh,ylow,yhigh,method1);
  nplot=nout;
  %[meqn,mx,my,method1,nplot,xl,yl,xc,yc,datafmt]=readgrid(outputdir);

  yflip=-1*[0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0];
  xflip=-1*[0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0];
  flip = -1*(yflip+2*xflip);
  %flip   =[0, 1, 2, 0, 0, 0, 1, 2, 0, 0, 2, 1, 3, 1, 2, 0, 3, 0];
 %flip_idx=[1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18];

  rho_i    =  1; %'rho_i';
  M1_i     =  2; %'M1_i';
  M2_i     =  3; %'M2_i';
  M3_i     =  4; %'M3_i';
  energy_i =  5; %'energy_i';
  rho_e    =  6; %'rho_e';
  M1_e     =  7; %'M1_e';
  M2_e     =  8; %'M2_e';
  M3_e     =  9; %'M3_e';
  energy_e = 10; %'energy_e';
  B1       = 11; %'B1';
  B2       = 12; %'B2';
  B3       = 13; %'B3';
  E1       = 14; %'E1';
  E2       = 15; %'E2';
  E3       = 16; %'E3';
  phi      = 17; %'phi';
  psi      = 18; %'psi';
  name = {'rho_i','M1_i','M2_i','M3_i','energy_i',...
          'rho_e','M1_e','M2_e','M3_e','energy_e',...
          'B1','B2','B3','E1','E2','E3','phi','psi'};

  % so that user can enter q to quit
  q=-1;
  qq=-2;
  
  % examples of valid responses:
  %  17
  %  phi
  %  [1,9,16]
  %  [rho_i,M3_e,E3]
  mf = 0;
  disp(' ')
  while(mf~=-1)
    % display choices
    for i=1:meqn
      disp([num2str(i), ' = ', char(name(i))]);
    end
    mf = input([ 'Which component[s] do you want to plot ( 1 - ',num2str(length(name)),...
            ', q=quit ) ? ']);
    if isempty(mf)
      mf=[1,9,16]; % default: density, current, E3
      %mf=[1:meqn];
    end
    if(mf<0)
      break;
    end

    for i=1:numel(mf)
      figure(mf(i));
      % leave the contents of the figure untouched
      % in case the user wants to use this feature
      % to make a comparison.
      %clf;
    end

    nf = 0;
    n1 = -1;
    while(nf~=-1)
      nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
            ' ) [type -1 or q to quit] ? ']);
      if isempty(nf)
        n1 = n1 + 1;
        nf = 0;
      else
        n1 = nf;
      end
      if n1> nplot
        disp(' ');
        disp(' End of plots ');
        disp(' ');
        n1 = nplot;
      end
      if(nf==qq)
        return;
      end
      if (nf==q)
        break;
      end
    
      [qfull,time]=readq2(n1,outputdir,datafmt,meqn,mx,my,method1);

      for i=1:length(mf)
         plot_scalar_outer(mf(i),qfull(:,:,mf(i)),xl,yl,time,...
            ['q(' num2str(mf(i)) ')=' char(name(mf(i)))],xflip(mf(i)),yflip(mf(i)));
      end
    end
  end
  disp(' ')
end

function plot_scalar_outer(fig,qvals,xl,yl,time,titleStr,xflip,yflip);

%    % bits of flip encode negation symmetry per axis
%    %
%    if(nargin<7)
%      flip=0;
%    end
%    xflip=1;
%    yflip=1;
%    if(flip==1)
%      xflip=-1;
%    elseif(flip==2)
%      yflip=-1;
%    elseif(flip==3)
%      xflip=-1;
%      yflip=-1;
%    end

    if(ishandle(fig))
      set(0,'CurrentFigure',fig);
    else
     figure(fig);
    end

    clf;

    color_range=[-0.2,0.2];
    plot_scalar(qvals,xl,yl,xflip,yflip,color_range);
%    qplot=zeros(size(qvals)+1);
%    qplot(1:end-1,1:end-1)=qvals;
%    pcolor(xl,yl,qplot);
%    shading flat;
%    colormap('jet');
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
%    shading flat;
    set_axes();
%    hold off;
%    axis on; box on; grid off;
%    axis('equal');
%    axis([-4*pi 4*pi -2*pi 2*pi]);
%    set(gca,'xtick',-12:4:12);
%    set(gca,'ytick',-6:2:6);
%    set(gca,'fontsize',16);
    set_title(titleStr,time);
%    set(t1,'fontsize',16);
%    c1=colorbar;
end
    
function plot_pseudovector(fig,xc,yc,qx,qy,time,titleStr);
    %qvec=qfull(:,:,11:12);
    %titleStr='B';
    %fig=40;
    if(ishandle(fig))
      set(0,'CurrentFigure',fig);
    else
      figure(fig);
    end
    clf;
    hold on;
    scale=.3;
    quiver(xc,yc,qx,qy,scale);
    quiver(xc,-yc,-qx,qy,scale);
    quiver(-xc,yc,qx,-qy,scale);
    quiver(-xc,-yc,-qx,-qy,scale);
    hold off;
    axis on; box on; grid off;
    axis('equal');
    axis([-4*pi 4*pi -2*pi 2*pi]);
    set(gca,'xtick',-12:4:12);
    set(gca,'ytick',-6:2:6);
    set(gca,'fontsize',16);
    t1 = title([titleStr ' at t = ',num2str(time)]); 
    set(t1,'fontsize',16);
end

