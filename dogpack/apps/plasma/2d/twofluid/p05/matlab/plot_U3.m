function plot_U3(outputdir_in,maxtime)

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

filename=[outputdir '/recon_E3andU3.dat'];
disp(['reading ' filename]);
fid = fopen(filename);
consv = fscanf(fid,'%e',[3 inf]);
status = fclose(fid);

t = consv(1,:)';
E3 = consv(2,:)';
U3 = consv(3,1)-consv(3,:)';
nt = length(t);


figure(101)
clf
hold on;
% until I fix my data display the reconnected flux 
% that the new algorithm would calculate by subtracting
% off the left_lost_flux
p1=plot(t,0.5*U3,'k-');
set(p1,'linewidth',2);
t1=title([ '0.5*(U3(0,0,0)-U3(0,0,t)) ' ...
         '(mx=' num2str(mx) ', my=' num2str(my) ')']);
set(t1,'fontsize',16);
axis on; box on; grid on;
hold off
if(nargin>=2)
   axis([0,maxtime,0,7]);
end
set(gca,'ytick',0:0.5:7.0);
set(gca,'plotboxaspectratio',[1.5 1 1]);

