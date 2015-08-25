function plot_resistivity(outputdir_in)
  global outputdir
  if(nargin<1)
    outputdir='output';
  else
    outputdir=outputdir_in;
  end
  
  format long e;
  
  % TITLE INFO --------------------------------------------------------------
  disp(' ')
  disp('  DOGPACK 2D GRID PLOTTING ROUTINE')
  disp(' ')
  % END TITLE INFO ----------------------------------------------------------
  
  %%% read in data
  
  % extract application parameters
  %
  %clear global mass_ratio;
  global mass_ratio;
  read_config_file([outputdir '/param.data']);

  % extract dogpack parameters
  global nout meqn mx my xlow xhigh ylow yhigh method datafmt;
  read_config_file([outputdir '/dogpack.data']);
  %[nout,tfinal,dtv,cflv,nv,...
  %   method,meqn,mx,my,mbc,xlow,xhigh,ylow,yhigh,...
  %   mrestart,nstart,datafmt]...
  % = read_dogpack_parameters2([outputdir '/dogpack.data']);
  [dx,dy,xl,yl,xc,yc]=get_grid(mx,my,xlow,xhigh,ylow,yhigh,method(1));
  nplot=nout;
  
  flip    =[0, 1, 2, 0, 0, 0, 1, 2, 0, 0, 2, 1, 2, 1, 2, 0, 3, 0];
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
    if(mf==-1)
      break;
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
      if (nf==-1)
        break
      end

      [qfull,time]=readq(n1,outputdir,datafmt,meqn,mx,my,method(1));
    
      %basename = sprintf('%04d',n1);
      %fname = [outputdir '/q',basename,'.dat'];
    
      %fids = fopen(fname,'r');
      %if(fids==-1)
      %    error(['could not open file ' fname]);
      %end
      %time = fscanf(fids,'%e',1);
      %qsoln = fscanf(fids,'%e',[meqn,inf]);
      %fclose(fids);
      %qsoln = transpose(qsoln);
      %
      %for me=1:meqn
      %  qfull(:,:,me) = reshape(qsoln(:,me),mx,my);
      %end
      mx_out = mx*method(1);
      my_out = my*method(1);
      rho_aug = zeros(mx_out+1,my_out+1);
      rho_aug(1:mx_out,1:my_out) = qfull(:,:,1) + qfull(:,:,6);
      
      % Construct fake resistivity
      u1tmp = (qfull(:,:,2)+qfull(:,:,7))./(qfull(:,:,1)+qfull(:,:,6));
      u2tmp = (qfull(:,:,3)+qfull(:,:,8))./(qfull(:,:,1)+qfull(:,:,6));
      u3tmp = (qfull(:,:,4)+qfull(:,:,9))./(qfull(:,:,1)+qfull(:,:,6));
      
      B1tmp = qfull(:,:,11);
      B2tmp = qfull(:,:,12);
      B3tmp = qfull(:,:,13);
      
      E1tmp = qfull(:,:,14);
      E2tmp = qfull(:,:,15);
      E3tmp = qfull(:,:,16);
      
      r = mass_ratio;
      
      J1tmp = (1+1/r)*qfull(:,:,2) - (1+r)*qfull(:,:,7);
      J2tmp = (1+1/r)*qfull(:,:,3) - (1+r)*qfull(:,:,8);
      J3tmp = (1+1/r)*qfull(:,:,4) - (1+r)*qfull(:,:,9);
      
      BxU1 = B2tmp.*u3tmp - B3tmp.*u2tmp;
      BxU2 = B3tmp.*u1tmp - B1tmp.*u3tmp;
      BxU3 = B1tmp.*u2tmp - B2tmp.*u1tmp;
      
      eta = (J1tmp.*(E1tmp-BxU1) + J2tmp.*(E2tmp-BxU2) + J3tmp.* ...
             (E3tmp-BxU3));
      eta = eta./(J1tmp.*J1tmp + J2tmp.*J2tmp + J3tmp.*J3tmp);

      Nonmag = sqrt((E1tmp-BxU1).^2 + (E2tmp-BxU2).^2 + (E3tmp- ...
                                                      BxU3).^2);
      Emag = sqrt(E1tmp.^2 + E2tmp.^2 + E3tmp.^2);
      BxUmag = sqrt(BxU1.^2 + BxU2.^2 + BxU3.^2);
      
      Jmag = sqrt(J1tmp.*J1tmp + J2tmp.*J2tmp + J3tmp.*J3tmp);
      
      eta_aug = zeros(mx_out+1,my_out+1);
      eta_aug(1:mx_out,1:my_out) = eta;
      
      Nonmag_aug = zeros(mx_out+1,my_out+1);
      Nonmag_aug(1:mx_out,1:my_out) = Nonmag;
      
      BxUmag_aug = zeros(mx_out+1,my_out+1);
      BxUmag_aug(1:mx_out,1:my_out) = BxUmag;
      
      Emag_aug = zeros(mx_out+1,my_out+1);
      Emag_aug(1:mx_out,1:my_out) = Emag;
      
      Jmag_aug = zeros(mx_out+1,my_out+1);
      Jmag_aug(1:mx_out,1:my_out) = Jmag;
      
      
      figure(1);
      clf;

      pcolor(xl,yl,log(rho_aug));
      hold on;
      pcolor(-xl,yl,log(rho_aug));
      pcolor(-xl,-yl,log(rho_aug));
      pcolor(xl,-yl,log(rho_aug));

      colormap('jet');
      shading flat;
      hold on;
      %p1=plot([-pi pi pi -pi -pi],[-pi/2 -pi/2 pi/2 pi/2 -pi/2],'k-');
      %set(p1,'linewidth',1);
      hold off;
      axis on; box on; grid off;
      axis('equal');
      axis([-4*pi 4*pi -2*pi 2*pi]);
      set(gca,'xtick',-12:4:12);
      set(gca,'ytick',-6:2:6);
      set(gca,'fontsize',16);
      
      t1 = title(['Log(Density) at t = ',num2str(time)]); 
      colorbar
      
      set(t1,'fontsize',16);
      
      
      figure(2);
      clf;

      pcolor(xl,yl,Nonmag_aug);
      hold on;
      pcolor(-xl,yl,Nonmag_aug);
      pcolor(-xl,-yl,Nonmag_aug);
      pcolor(xl,-yl,Nonmag_aug);

      colormap('jet');
      shading flat;
      hold on;
      %p1=plot([-pi pi pi -pi -pi],[-pi/2 -pi/2 pi/2 pi/2 -pi/2],'k-');
      %set(p1,'linewidth',1);
      hold off;
      axis on; box on; grid off;
      axis('equal');
      axis([-4*pi 4*pi -2*pi 2*pi]);
      set(gca,'xtick',-12:4:12);
      set(gca,'ytick',-6:2:6);
      set(gca,'fontsize',16);
      
      t1 = title(['|| E - BxU ||  at t = ',num2str(time)]); 
      colorbar;

      set(t1,'fontsize',16);
      
      
      figure(3);
      clf;
      
      pcolor(xl,yl,Emag_aug);
      hold on;
      pcolor(-xl,yl,Emag_aug);
      pcolor(-xl,-yl,Emag_aug);
      pcolor(xl,-yl,Emag_aug);

      colormap('jet');
      shading flat;
      hold on;
      %p1=plot([-pi pi pi -pi -pi],[-pi/2 -pi/2 pi/2 pi/2 -pi/2],'k-');
      %set(p1,'linewidth',1);
      hold off;
      axis on; box on; grid off;
      axis('equal');
      axis([-4*pi 4*pi -2*pi 2*pi]);
      set(gca,'xtick',-12:4:12);
      set(gca,'ytick',-6:2:6);
      set(gca,'fontsize',16);
      
      t1 = title(['|| E ||  at t = ',num2str(time)]); 
      colorbar;

      set(t1,'fontsize',16);
      
      
      figure(4);
      clf;
      
      pcolor(xl,yl,BxUmag_aug);
      hold on;
      pcolor(-xl,yl,BxUmag_aug);
      pcolor(-xl,-yl,BxUmag_aug);
      pcolor(xl,-yl,BxUmag_aug);

      colormap('jet');
      shading flat;
      hold on;
      %p1=plot([-pi pi pi -pi -pi],[-pi/2 -pi/2 pi/2 pi/2 -pi/2],'k-');
      %set(p1,'linewidth',1);
      hold off;
      axis on; box on; grid off;
      axis('equal');
      axis([-4*pi 4*pi -2*pi 2*pi]);
      set(gca,'xtick',-12:4:12);
      set(gca,'ytick',-6:2:6);
      set(gca,'fontsize',16);
      
      t1 = title(['|| BxU ||  at t = ',num2str(time)]); 
      colorbar;

      set(t1,'fontsize',16);
      
      
      figure(5);
      clf;
      
      pcolor(xl,yl,Jmag_aug);
      hold on;
      pcolor(-xl,yl,Jmag_aug);
      pcolor(-xl,-yl,Jmag_aug);
      pcolor(xl,-yl,Jmag_aug);

      colormap('jet');
      shading flat;
      hold on;
      %p1=plot([-pi pi pi -pi -pi],[-pi/2 -pi/2 pi/2 pi/2 -pi/2],'k-');
      %set(p1,'linewidth',1);
      hold off;
      axis on; box on; grid off;
      axis('equal');
      axis([-4*pi 4*pi -2*pi 2*pi]);
      set(gca,'xtick',-12:4:12);
      set(gca,'ytick',-6:2:6);
      set(gca,'fontsize',16);
      
      t1 = title(['|| J ||  at t = ',num2str(time)]); 
      colorbar;

      set(t1,'fontsize',16);
      
      
      
      % for elc-pos case
      %for i=1:length(mf)
      %   plot_scalar(mf(i),qfull(:,:,mf(i)),xl,yl,time,['q(' num2str(mf(i)) ')=' char(name(mf(i)))],flip(mf(i)));
      %end
    end
  end
  disp(' ')
end

function plot_scalar(fig,qvals,xl,yl,time,titleStr,flip);

    % bits of flip encode negation symmetry per axis
    %
    if(nargin<7)
      flip=0;
    end
    xflip=1;
    yflip=1;
    if(flip==1)
      xflip=-1;
    elseif(flip==2)
      yflip=-1;
    elseif(flip==3)
      xflip=-1;
      yflip=-1;
    end

    %if(ishandle(fig))
    %  set(0,'CurrentFigure',fig);
    %else
      figure(fig);
    %end

    clf;
    qplot=zeros(size(qvals)+1);
    qplot(1:end-1,1:end-1)=qvals;
    pcolor(xl,yl,qplot);
    shading flat;
    colormap('jet');
    hold on;
    if(yflip==-1) pcolor( xl,-yl,-qplot);
    else pcolor( xl,-yl,qplot);
    end;
    shading flat;
    if(xflip==-1) pcolor(-xl, yl,-qplot);
    else pcolor(-xl, yl,qplot);
    end
    shading flat;
    pcolor(-xl,-yl,(xflip*yflip)*qplot);
    shading flat;
    hold off;
    axis on; box on; grid off;
    axis('equal');
    axis([-4*pi 4*pi -2*pi 2*pi]);
    set(gca,'xtick',-12:4:12);
    set(gca,'ytick',-6:2:6);
    set(gca,'fontsize',16);
    t1 = title([titleStr ' at t = ',num2str(time) ' ']); 
    set(t1,'fontsize',16);
    c1=colorbar;
    set(c1,'fontsize',16);
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

function plot_in_current_pane(x,pane_data,low_bound,high_bound,...
     fnt,pane_name)

  global plot_lower plot_upper tick_increment output_cells_per_plot_cell

  plot_handle=plot(x,pane_data,'k-');
  hold on;
  set(gca,'fontsize',fnt);
  %set(gca,'fontweight','bold');
  set(gca,'linewidth',1);
  set(plot_handle,'linewidth',1);
  box on; grid off; axis on;
%  set(gca,'xtick',-0.4:0.2:0.4);
%  set(gca,'xlim',[-0.5 0.5]);
  %set(gca,'xtick',-2.0:0.2:0.0);
  %set(gca,'xlim',[-2.0 0.0]);

  % use data from setprob.data to define axes.
  set(gca,'xtick',plot_lower:tick_increment:plot_upper);
  set(gca,'xlim',[plot_lower plot_upper]);

  set(gca,'ylim',[low_bound high_bound]);
  t1 = title(pane_name);
  set(t1,'fontsize',fnt);
  %set(t1,'fontweight','bold');
  %ylabel(pane_name);
  set(gca,'fontsize',fnt);
  hold off;

end
