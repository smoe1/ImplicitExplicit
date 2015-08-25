%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%                   meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc): grid points (cell centers), size = (mx,my)
%%       (xl,yl): grid points (lower left cell corners), size = (mx+1,my+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  figure(1);
  clf;
  surf(xl,yl,qaug(:,:,m));
  colormap('jet');
  axis on; box on; grid off;
  axis('equal');
  axis auto;
  axis([xl(1,1) xl(mx,1) yl(1,1) yl(1,my)]);
  set(gca,'fontsize',16);
  t1 = title(['f(t,x,v) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);

  shading flat;
  colorbar;
  caxis auto;
  %caxis([0,0.6]);
  set(gca,'fontsize',16);
  set(gca,'xtick',-10:5:10);
  set(gca,'ytick',-6:2:6);


  %disp([['  minimum value of f is = ', num2str(min(min(qsoln)), '%2.8e')]])

  figure(2)
  clf;
  if mod(mx,2) == 0 % even number of cells printed!
    xcenter_index = mx / 2;
  else
    xcenter_index = (mx+1) / 2;
  end

  vel = yc(xcenter_index,:);
  qslice = qsoln(xcenter_index,:)';
  p1=plot(vel, qslice,'r-');
  set(p1,'linewidth',2);
  axis([yl(1,1) yl(1,my) -0.05 0.42]);
  hold on;
  p2=plot(vel, zeros(size(vel)), 'b-');  % show y = 0 on the grid
  set(p2,'linewidth',2);
  set(gca,'fontsize',16);
  %t1 = title(['f(x = ',num2str(xc(xcenter_index,5)),',v,t) at t = ',num2str(time),'     [DoGPack]']); 
  t1 = title(['f(t,x = ',num2str(xc(xcenter_index,5), '%1.3f'), ',v) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);

  figure(1)

  %%%% Print the pictures!
  picname1 = strcat(outputdir, '_t', num2str(time), '.eps'  );
  picname2 = strcat(outputdir, '_t', num2str(time), '_slice', '.eps');
% 
  print(1, '-depsc2', picname1);
  print(2, '-depsc2', picname2);
  
  
  % OUTPUT TO FILE
% if abs(time)<1.0e-14
%   ind = 0;
%   fidgrid = fopen(['bump_on_tail/mesh_info.dat'],'wt');
%   fprintf(fidgrid,'%10d\n',mx);
%   fprintf(fidgrid,'%10d\n',my);
%   fprintf(fidgrid,'%24.16e\n',xlow);
%   fprintf(fidgrid,'%24.16e\n',xhigh);
%   fprintf(fidgrid,'%24.16e\n',ylow);
%   fprintf(fidgrid,'%24.16e\n',yhigh);
%   fprintf(fidgrid,'%10d\n',nplot);
%   fclose(fidgrid);
% else
%   ind = ind+1;
% end
% 
% if (ind<10)
%   fidout = fopen(['bump_on_tail/f0',num2str(ind),'.dat'],'wt');
% else
%   fidout = fopen(['bump_on_tail/f',num2str(ind),'.dat'],'wt');
% end
% fprintf(fidout,'%24.16e\n',time);
% for j=1:my
%   for i=1:mx
%     fprintf(fidout,'%24.16e\n',qsoln(i,j,1));    
%   end
% end
% fclose(fidout);
