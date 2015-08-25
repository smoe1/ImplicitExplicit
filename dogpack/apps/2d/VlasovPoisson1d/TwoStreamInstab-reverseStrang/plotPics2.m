function plotPics2( outputdir )

  mcons = 5;  % number of conserved quantities 

  for n=1:5
    figure(n);
    clf;
  end

  fname = [outputdir '/Efield.dat'];
  
  fid = fopen(fname);
  if fid==-1
    error(['File  ',outputdir,'  not found.']);
  end
  consv = fscanf(fid,'%e',[mcons+1 inf]);

  status = fclose(fid);

  t  = consv(1,:)';
  nt = length(t);
  tmax = t(end);

  % Electric Field
  normE = consv(2,:)';

  figure(1);

  q = semilogy(t,normE,'ro');
  hold on;

  % L1 norm of f
  figure(2);

  qc   = consv(3,:)';
  dqc  = (qc - qc(1))/qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % L2 norm of f
  figure(3);
  qc   = consv(4,:)';
  dqc  = (qc - qc(1))/qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % Total Energy
  figure(4);
  qc   = consv(5,:)';
  dqc  =  (qc - qc(1))/qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % Entropy 
  figure(5);
  qc   = consv(6,:)';
  dqc  =  (qc - qc(1))/qc(1);
  p = plot(t, dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % add titles:
  figure(1);
  set(gca,'fontsize',16);
  t1=title([ 'L_2(Electric Field)']);
  set(t1,'fontsize',16);
  axis([0 tmax 1e-2 1e1]);
  yticks = [1e-2 1e-1 1e0 1e1];
  set(gca, 'ytick', yticks );
  %l = legend('SL2', 'SL4', 'location', 'north');

  % L1 norm
  figure(2);
  set(gca,'fontsize',16);
  t1=title([ '(L_1( f(t) ) - L_1( f(t=0) ))/L_1( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 tmax -1.5e-14 1.5e-14]);
  set(gca,'xtick',0:5:45);
  set(gca,'ytick',-1.5e-14:0.5e-14:1.5e-14);
  %l = legend('SL2', 'SL4', 'location', 'north');

  % L2 norm
  figure(3);
  set(gca,'fontsize',16);
  t1=title([ '(L_2( f(t) ) - L_2( f(t=0) ))/L_2( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.07 0.005]);
  set(gca,'xtick',0:5:45);
  set(gca, 'ytick', [-0.07:0.01:0]);
  %l = legend('SL2', 'SL4', 'location', 'southwest');

  % energy
  figure(4);
  set(gca,'fontsize',16);
  t1=title(['(Energy( t ) - Energy( t=0 ))/Energy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.0005 0.00125]);
  set(gca,'xtick',0:5:45);
  set(gca,'ytick', [-0.0005:0.00025:0.00125]);
  %l = legend('SL2', 'SL4', 'location', 'northwest');

  % entropy
  figure(5);
  set(gca,'fontsize',16);
  t1=title(['(Entropy( t ) - Entropy( t=0 ))/Entropy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.01 0.04]);
  set(gca,'xtick',0:5:45);
  set(gca, 'ytick', [0.0:0.01:0.04]);
  %l = legend('SL2', 'SL4', 'location', 'northwest');
  %set(l, 'fontsize', 16 );

  for n=1:5
    figure(n);
    set(gca, 'xtick', 0:5:tmax);
  end

end
