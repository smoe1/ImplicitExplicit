function plot_recon(outputdir_in,maxtime,fig)

  global outputdir;
  if(nargin>=1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  if(~isdir(outputdir))
    error(['no such directory: ' outputdir]);
  end
  if(nargin<3)
     fig=102;
  if(nargin<2)
     maxtime=40;
  end
  end

format long e

% INFO --------------------------------------------------------------
disp(' ')
disp('Reconnection Plotter')
disp(' ')

% INFO --------------------------------------------------------------

% get information for plotting purposes
global mx my mass_ratio temp_ratio elc_mass ion_mass;
read_config

disp(['plotting ' outputdir ...
      ', elc_mass=' num2str(elc_mass) ...
      ', ion_mass=' num2str(ion_mass) ]);

filename1=[outputdir '/rec_flx.dat'];
filename2=[outputdir '/Xpoint.dat'];
%filename=[outputdir '/recon_flux.dat'];
disp(['reading ' filename1 ' and ' filename2]);
fid= fopen(filename1);
consv1 = fscanf(fid,'%e',[6 inf]);
status = fclose(fid);
fid= fopen(filename2);
consv2 = fscanf(fid,'%e',[7 inf]);
status = fclose(fid);

t = consv1(1,:)';
recon_flux     = consv1(2,:)';
rght_flux      = consv1(3,:)';
left_flux      = consv1(4,:)';
left_lost_flux = consv1(5,:)';
total_x_axis_island_flux = consv1(6,:)';
 E3           = consv2(2,:)';
 E3_accum     = consv2(3,:)';
 J3           = consv2(4,:)';
 rho          = consv2(5,:)';
 J3_ovr_rho   = consv2(6,:)'; %J3./rho;
 left_flux    = consv2(7,:)';
nt = length(t);

mime_J3_ovr_rho = J3./rho*(ion_mass*elc_mass);

figure(fig)
clf;
hold on;
% until I fix my data display the reconnected flux 
% that the new algorithm would calculate by subtracting
% off the left_lost_flux
%p1=plot(t,recon_flux-left_lost_flux,'k-');
p1=plot(t,left_flux(1)-left_flux,'r-');
set(p1,'linewidth',3);
p1=plot(t,-E3_accum,'b--');
set(p1,'linewidth',3);
p1=plot(t,-mime_J3_ovr_rho+mime_J3_ovr_rho(1),'k-.');
set(p1,'linewidth',3);
p1=plot(t,rght_flux,'b-.');
set(p1,'linewidth',2);
p1=plot(t,left_flux,'g-'); 
set(p1,'linewidth',2);
p1=plot(t,left_lost_flux,'r--'); % size of central island
set(p1,'linewidth',2);
p1=plot(t,total_x_axis_island_flux-left_lost_flux,'m-.'); % size of other islands
set(p1,'linewidth',2);
p1=plot(t,rho,'m-');
set(p1,'linewidth',2);
hold off;
set(gca,'fontsize',16);
%t1=title([ 'X-point values and boundary fluxes' ]);
t1=title([ 'Reconnection for' ...
         ' m_i/m_e=' num2str(mass_ratio) ...
         ', T_i/T_e=' num2str(temp_ratio) ...
         ' (mx=' num2str(mx) ', my=' num2str(my) ...
         ')']);
set(t1,'fontsize',32);
axis on; box on; grid on;
hold off
axis([0,maxtime,-0.5,7]);
set(gca,'xtick',0:5:50);
set(gca,'ytick',0:0.5:7.0);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'YAxisLocation','right');
set(gca,'fontsize',24);
l1=legend(...
   'reconnected flux', ...
   'time integral of (-E3)', ...
   '\Delta (-mi*me)J3/rho', ...
   'flux exiting right wall', ...
   'flux entering left wall', ...
   'flux of central island', ...
   'flux of other islands', ...
   'mass density' );
%set(l1,'location','northeast');
set(l1,'location','west');
set(l1,'fontsize',24);

