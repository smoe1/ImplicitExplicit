function plot_E( mxsize )

  mcons = 5;  % number of conserved quantities 

  for n=1:5
    figure(n);
    clf;
  end

  fid = fopen([strcat('outputSL2_', int2str(mxsize)) '/Efield.dat']);
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
  dqc  = qc - qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % L2 norm of f
  figure(3);
  qc   = consv(4,:)';
  dqc  = qc - qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % Total Energy
  figure(4);
  qc   = consv(5,:)';
  dqc  = qc - qc(1);
  p = plot(t,dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  % Entropy 
  figure(5);
  qc   = consv(6,:)';
  dqc  = qc - qc(1);
  p = plot(t, dqc,'r--');
  set(p, 'linewidth', 3 );
  hold off;

  %%%%%%%%%%%% Second plot on same figures %%%%%%%%%%%%%%%%%%%%
  fid = fopen([strcat('outputSL4_', int2str(mxsize)) '/Efield.dat']);
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
  q = semilogy(t,normE,'b-', 'linewidth', 3);
  hold off;

  % L1 norm of f
  figure(2);
  qc   = consv(3,:)';
  dqc  = qc - qc(1);
  hold on;
  p = plot(t,dqc,'b-', 'linewidth', 3);
  hold on;
  q = plot(t, zeros(size(t)), 'k-', 'linewidth', 3);
  hold off;

  % L2 norm of f
  figure(3);
  qc   = consv(4,:)';
  dqc  = qc - qc(1);
  hold on;
  p = plot(t,dqc,'b-');
  set(p, 'linewidth', 3 );
  hold on;
  q = plot(t, zeros(size(t)), 'k-' );
  set(q, 'linewidth', 3 );
  hold off;

  % Total Energy
  figure(4);
  qc   = consv(5,:)';
  dqc  = qc - qc(1);
  hold on;
  p = plot(t,dqc,'b-', 'linewidth', 3 );
  hold on;
  q = plot(t, zeros(size(t)), 'k-' );
  set(q, 'linewidth', 3 );
  hold off;

  % Entropy 
  figure(5);
  qc   = consv(6,:)';
  dqc  = qc - qc(1);
  hold on;
  p = plot(t, dqc,'b-');
  set(p, 'linewidth', 3 );
  hold on;
  q = plot(t, zeros(size(t)), 'k-' );
  set(q, 'linewidth', 3 );
  hold off;

  % add titles:
  figure(1);
  set(gca,'fontsize',16);
  t1=title([ 'L_2(Electric Field)']);
  set(t1,'fontsize',16);
  axis([0 tmax 1e-2 1e1]);
  yticks = [1e-2 1e-1 1e0 1e1];
  set(gca, 'ytick', yticks );
  l = legend('SL2', 'SL4', 'location', 'north');

  % L1 norm
  figure(2);
  set(gca,'fontsize',16);
  t1=title([ 'L_1( f(t) ) - L_1( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 tmax -1e-12 1e-12]);
  l = legend('SL2', 'SL4', 'location', 'north');

  % L2 norm
  figure(3);
  set(gca,'fontsize',16);
  t1=title(['L_2( f(t) ) - L_2( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.2 0.05]);
  set(gca, 'ytick', [-0.2:0.05:0.1]);
  l = legend('SL2', 'SL4', 'location', 'southwest');

  % energy
  figure(4);
  set(gca,'fontsize',16);
  t1=title(['Energy( t ) - Energy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.01 0.035]);
  set(gca, 'ytick', [-0.01:0.005:.035]);
  l = legend('SL2', 'SL4', 'location', 'northwest');

  % entropy
  figure(5);
  set(gca,'fontsize',16);
  t1=title(['Entropy( t ) - Entropy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 tmax -0.1 1.0]);
  set(gca, 'ytick', [-0.1:0.1:1]);
  l = legend('SL2', 'SL4', 'location', 'northwest');
  set(l, 'fontsize', 16 );

  for n=1:5
    figure(n);
    set(gca, 'xtick', 0:5:tmax);
  end

  % print all the pretty pictures
  name1 = strcat( 'Two_', strcat( int2str(mxsize), '_ElectricField.eps' ) );
  name2 = strcat( 'Two_', strcat( int2str(mxsize), '_L1f.eps' ) );
  name3 = strcat( 'Two_', strcat( int2str(mxsize), '_L2f.eps' ) );
  name4 = strcat( 'Two_', strcat( int2str(mxsize), '_Energy.eps' ) );
  name5 = strcat( 'Two_', strcat( int2str(mxsize), '_Entropy.eps' ) );
  print(1, '-depsc2', name1 );
  print(2, '-depsc2', name2 );
  print(3, '-depsc2', name3 );
  print(4, '-depsc2', name4 );
  print(5, '-depsc2', name5 );

end
