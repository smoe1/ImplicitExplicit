%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get MHD variables
[xc_mhd,qsoln_mhd]=get_mhd(n1,points_per_dir,'mhd_compare/output');
p_mhd =  (5/3-1)*(qsoln_mhd(:,5) ...
                  -0.5*(qsoln_mhd(:,2).^2+qsoln_mhd(:,3).^2+ qsoln_mhd(:,4).^2)./qsoln_mhd(:,1)...
                  -0.5*(qsoln_mhd(:,6).^2+qsoln_mhd(:,7).^2+ qsoln_mhd(:,8).^2));

% parameters
fids = fopen('param.data','r');
gamma      = fscanf(fids,'%g',1);   fscanf(fids,'%s',1); ...
    fscanf(fids,'%s',1); fscanf(fids,'%s',1);
mass_ratio = fscanf(fids,'%g',1);   fscanf(fids,'%s',1); fscanf(fids,'%s',1);
debye      = fscanf(fids,'%g',1);   fscanf(fids,'%s',1); fscanf(fids,'%s',1);
cs_light   = fscanf(fids,'%g',1);   fscanf(fids,'%s',1); ...
    fscanf(fids,'%s',1); fscanf(fids,'%s',1);
larmor     = fscanf(fids,'%g',1);   fscanf(fids,'%s',1);
fclose(fids);

figure(1);
clf;
pm=plot(xc_mhd,qsoln_mhd(:,1),'r-');
hold on;
pz=plot(xc,qsoln(:,1)+qsoln(:,6),'b-');
set(pm,'linewidth',1);
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.2:1.2);
axis([0 1 0 1.2]);
t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']);
    
        
press_i = (5/3-1)*(qsoln(:,5) -0.5*(qsoln(:,2).^2+qsoln(:,3).^2+ ...
                                    qsoln(:,4).^2)./qsoln(:,1));

press_e = (5/3-1)*(qsoln(:,10)-0.5*(qsoln(:,7).^2+qsoln(:,8).^2+ ...
                                    qsoln(:,9).^2)./qsoln(:,6));
    
    
figure(2);
clf;
pm=plot(xc_mhd,p_mhd,'r-');
hold on;
pz=plot(xc,press_i + press_e,'b-');
set(pm,'linewidth',1);
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.2:1.2);
axis([0 1 0 1.3]);
t1 = title(['Pressure at t = ',num2str(time),'     [DoGPack]']);
      
figure(3);
clf;
pm=plot(xc_mhd,qsoln_mhd(:,7),'r-');
hold on;
pz=plot(xc,qsoln(:,12),'b-');
set(pm,'linewidth',1);
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.5:0.5:1.5);
axis([0 1 -1.5 1.5]);
t1 = title(['B^2 at t = ',num2str(time),'     [DoGPack]']);
