clear;

fids = fopen('shllw.dat','r');
tmp = fscanf(fids,'%e %e %e %e',[4,inf]);
fclose(fids);
x = transpose(tmp(1,:));
h = transpose(tmp(2,:));
b = transpose(tmp(3,:));
Fr = transpose(tmp(4,:));
clear tmp;

figure(1)
clf
plot(x,b+h,'k-');
hold on;
fill(x,b,'g');
hold off;
set(gca,'xlim',[0 1]);
set(gca,'ylim',[0 1.1*max(b+h)]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1.2*max(b+h));
set(gca,'fontsize',16);
t1=title(['Steady State Solution for Shallow Water']);
