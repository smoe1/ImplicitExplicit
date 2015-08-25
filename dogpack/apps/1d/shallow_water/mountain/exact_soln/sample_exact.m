clear;

disp(' ')
display(' Loading grid data: please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% READ-IN GEOMETRIC INFO %%%

fids = fopen('../grid/elements2.dat','r');
tmp = fscanf(fids,'%d %d %d %e',[4,inf]);
fclose(fids);
t = transpose(tmp(1:3,:));
Area = transpose(tmp(4,:));
clear tmp;

fids = fopen('../grid/nodes2.dat','r');
tmp = fscanf(fids,'%e',[3,inf]);
fclose(fids);
p = transpose(tmp(1:2,:));
Cdual = transpose(tmp(3,:));
clear tmp;

fids = fopen('../grid/boundary2.dat','r');
tmp = fscanf(fids,'%e',[1,inf]);
fclose(fids);
bnd = transpose(tmp(1:1,:));
clear tmp;
for i=1:length(bnd)
  bnd_node(i,1:2) = p(bnd(i),1:2);
end 

display(' Finished loading grid data.');
disp(' ')


fids = fopen('shllw.dat','r');
tmp = fscanf(fids,'%e %e %e %e %e %e',[6,inf]);
fclose(fids);
xex(:,1) = transpose(tmp(1,:));
qex(:,1) = transpose(tmp(2,:));
qex(:,2) = transpose(tmp(3,:));
qex(:,3) = transpose(tmp(4,:));
bot_ex(:,1) = transpose(tmp(5,:));
Fr_ex(:,1)  = transpose(tmp(6,:));
clear tmp;

qint(:,1) = interp1(xex(:,1),qex(:,1),p(:,1),'pchip','extrap');
qint(:,2) = interp1(xex(:,1),qex(:,2),p(:,1),'pchip','extrap');
qint(:,3) = interp1(xex(:,1),qex(:,3),p(:,1),'pchip','extrap');

fids_ex = fopen('interp_ex.dat','w');
for i=1:length(p)
  fprintf(fids_ex,'%27.20e %27.20e %27.20e %27.20e\n',[p(i,1),qint(i,1),qint(i,2),qint(i,3)]);
end
fclose(fids_ex);


