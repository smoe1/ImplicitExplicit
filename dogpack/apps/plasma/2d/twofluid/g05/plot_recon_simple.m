function plot_recon_simple(outputdir_in,maxtime)

  global outputdir;
  if(nargin>=1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  if(~isdir(outputdir))
    error(['no such directory: ' outputdir]);
  end

format long e

% INFO --------------------------------------------------------------
disp(' ')
disp('Reconnection Plotter')
disp(' ')

% INFO --------------------------------------------------------------

% get information for plotting purposes
global mx my mass_ratio temp_ratio;
read_config

%filename=[outputdir '/rec_flx.dat'];
filename=[outputdir '/recon_flux.dat'];
disp(['reading ' filename]);
fid = fopen(filename);
consv = fscanf(fid,'%e',[6 inf]);
status = fclose(fid);

t = consv(1,:)';
recon_flux = consv(2,:)';
rght_flux = consv(3,:)';
left_flux = consv(4,:)';
left_lost_flux = consv(5,:)';
total_x_axis_island_flux = consv(6,:)';
nt = length(t);


figure(100)
clf
hold on;
% until I fix my data display the reconnected flux 
% that the new algorithm would calculate by subtracting
% off the left_lost_flux
p1=plot(t,recon_flux-left_lost_flux,'k-');
set(p1,'linewidth',2);
%p1=plot(t,rght_flux,'b-');
%p1=plot(t,left_flux,'g-'); 
%p1=plot(t,left_lost_flux,'r-'); % size of central island
%p1=plot(t,total_x_axis_island_flux-left_lost_flux,'m-'); % size of other islands
%flux_check_rght_flux = left_flux-left_lost_flux + recon_flux;
%p1=plot(t,flux_check_rght_flux,'y:'); % should agree with rght_flux
hold off;
set(gca,'fontsize',16);
t1=title([ 'Reconnected flux vs. time ' ]);% ...
%         '(mx=' num2str(mx) ', my=' num2str(my) ...
%         ', m_i/m_e=' num2str(mass_ratio) ...
%         ', T_i/T_e=' num2str(temp_ratio) ')']);
set(t1,'fontsize',16);
axis on; box on; grid on;
hold off
if(nargin>=2)
   axis([0,maxtime,0,7]);
   %set(gca,'xtick',0:100:500);
end
set(gca,'ytick',0:0.5:7.0);
set(gca,'plotboxaspectratio',[1.5 1 1]);

