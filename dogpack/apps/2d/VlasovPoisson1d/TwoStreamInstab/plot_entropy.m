clear
clear functions

format long e

outputdir = 'output';

% Parse qhelp.dat
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=2)
    error(['Incorrect dimension, ndims must be 2. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridTmp = GridType(1:9);
if (GridTmp=='Cartesian')
  GridType='Cartesian   ';
elseif (GridTmp=='Unstructu')
  GridType='Unstructured';
else
  error(['Incorrect GridType, must be Cartesian or Unstructured. GridType = ',GridType]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CARTESIAN plotting routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (GridType=='Cartesian   ')
    
  meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  meth1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  datafmt = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  my      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  ylow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  yhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
end
fclose(fids);

m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
	' ) ? ']);
disp(' ')
if isempty(m)
  m=1;
end

% INFO --------------------------------------------------------------
fid    = fopen('output/Efield.dat');
consv  = fscanf(fid,'%e',[6,inf]  );
t      = consv(1,:)';
nt     = length(t);

% L2-norm of electric field, l1-norm of solution, l2-norm of solution, Energy and Entroy
E2   = consv(2,:)';      E2diff  = E2 - E2(1);
l1   = consv(3,:)';      l1diff  = l1 - l1(1);
l2   = consv(4,:)';      l2diff  = l2 - l2(1);
En   = consv(5,:)';      Endiff  = En - En(1);
Ent  = consv(6,:)';      Entdiff = Ent - Ent(1);
status = fclose(fid);


figure(1);
clf;
plot( t, Entdiff );
set(gca,'fontsize',16);
t1 = title([ 'Conservation of Entropy vs. time']);
set(t1,'fontsize',16);
axis on; box on; grid off;
hold off

