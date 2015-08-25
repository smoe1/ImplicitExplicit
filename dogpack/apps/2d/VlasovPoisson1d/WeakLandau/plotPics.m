function plot_E( mxsize )

  mcons = 5;  % number of conserved quantities 

  fid = fopen([strcat('outputSL2_', int2str(mxsize)) '/Efield.dat']);
% fid = fopen([strcat('outputRK4_', int2str(mxsize)) '/Efield.dat']);
  if fid==-1
    error(['File  ',outputdir,'  not found.']);
  end
  consv = fscanf(fid,'%e',[mcons+1 inf]);

  status = fclose(fid);

  t  = consv(1,:)';
  nt = length(t);

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

  fid = fopen([strcat('outputSL4_', int2str(mxsize)) '/Efield.dat']);
  if fid==-1
    error(['File  ',outputdir,'  not found.']);
  end
  consv = fscanf(fid,'%e',[mcons+1 inf]);
  status = fclose(fid);

  t  = consv(1,:)';
  nt = length(t);


  % Electric Field
  normE = consv(2,:)';

  figure(1);
  q = semilogy(t,normE,'b-', 'linewidth', 3);
  hold on;

  %% add in lines lines for plot %
  gamma = -0.1533;
 
  if mxsize == 64
    y  = consv(1,2) * exp( gamma*(t) );
  else
    y  = (consv(1,2)+1e-2) * exp( gamma*(t) );
  end
  y = 0.06 * exp( gamma*t );

  p = semilogy(t,y,'k-', 'linewidth', 3);

  hold off;

  % L1 norm of f
  figure(2);

  %axis([0 60 -1e-12 1e-12]);
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
% p = plot(t,dqc,'b-', 'linewidth', 3 );
  p = plot(t,dqc,'bo', 'linewidth', 3 );
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
  axis([0 60 1e-10 1e-1]);

  yticks = 10.^( -(10:-1:1) );
  set(gca, 'ytick', yticks);
  set(gca, 'xtick', 0:5:60 );

  l = legend('SL2', 'SL4' );
  set(l, 'location', 'northeast');

  % L1 norm
  figure(2);
  set(gca,'fontsize',16);
  t1=title([ 'L_1( f(t) ) - L_1( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 60 -1e-12 1e-12]);
  l = legend('SL2', 'SL4', 'location', 'north');

  % L2 norm
  figure(3);
  set(gca,'fontsize',16);
  t1=title(['L_2( f(t) ) - L_2( f(t=0) )']);
  set(t1,'fontsize',16);
  axis([0 60 -15e-9 5e-9]);
  l = legend('SL2', 'SL4', 'location', 'southwest');

  % energy
  figure(4);
  set(gca,'fontsize',16);
  %t1=title(['Energy']);
  t1=title(['Energy( t ) - Energy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 60 -5e-7 15e-7]);
  l = legend('SL2', 'SL4', 'location', 'northeast');

  % entropy
  figure(5);
  set(gca,'fontsize',16);
  %t1=title(['Entropy']);
  t1=title(['Entropy( t ) - Entropy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 60 -2e-8 6e-8]);
  l = legend('SL2', 'SL4', 'location', 'northwest');
  set(l, 'fontsize', 16 );


  % print all the pretty pictures
% name1 = strcat( 'Weak_', strcat( int2str(mxsize), '_ElectricField.eps' ) );
% name2 = strcat( 'Weak_', strcat( int2str(mxsize), '_L1f.eps' ) );
% name3 = strcat( 'Weak_', strcat( int2str(mxsize), '_L2f.eps' ) );
% name4 = strcat( 'Weak_', strcat( int2str(mxsize), '_Energy.eps' ) );
% name5 = strcat( 'Weak_', strcat( int2str(mxsize), '_Entropy.eps' ) );
% print(1, '-depsc2', name1 );
% print(2, '-depsc2', name2 );
% print(3, '-depsc2', name3 );
% print(4, '-depsc2', name4 );
% print(5, '-depsc2', name5 );

end
