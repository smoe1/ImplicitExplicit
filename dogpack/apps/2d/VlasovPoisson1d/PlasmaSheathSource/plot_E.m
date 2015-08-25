function plot_E(outputdir_in)

  global outputdir
  if(nargin<1)
    outputdir='output';
  else
    outputdir=outputdir_in;
  end

  format long e

  % meqn and nplot
  fids  = fopen([outputdir,'/qhelp.dat'],'r');
  if fids==-1
    error(['File  ',outputdir,'/qhelp.dat  not found.']);
  end

  %%%%%%%%%%%%%% Read in all the parameters %%%%%%%%%%%%%
  JUNK  = fscanf(fids,'%s',1);
  meqn  = fscanf(fids,'%d',1);
  maux  = fscanf(fids,'%d',1);
  nplot = fscanf(fids,'%d',1);
  meth1 = fscanf(fids,'%d',1);
  mx    = fscanf(fids,'%d',1);
  my    = fscanf(fids,'%d',1);
  xlow  = fscanf(fids,'%e',1);
  xhigh = fscanf(fids,'%e',1);
  ylow  = fscanf(fids,'%e',1);
  yhigh = fscanf(fids,'%e',1);
  datafmt = fscanf(fids,'%e',1);
  fclose(fids);
 
  m = 1;

  % clear all the current figures:
  for n=1:5
    figure(n);
    clf
  end

  % INFO -------------------------------------------------------------- %

  mcons = 5;  % number of conserved quantities 

  %fid = fopen([outputdir '/Efield.dat']);
% fid = fopen(['outputSL2_64/Efield.dat']);
  fid = fopen(['outputSL2_128/Efield.dat']);
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
  set(gca,'fontsize',16);
  t1=title([ 'Strong Landau Damping: Electric Field vs. Time']);
  set(t1,'fontsize',16);

  q = semilogy(t,normE,'ro');
  hold on;

  % L1 norm of f
  figure(2);

  axis on; box on; grid off;
  %axis([0 60 -0.001 0.001]);
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

% fid = fopen(['outputSL4_64/Efield.dat']);
  fid = fopen(['outputSL4_128/Efield.dat']);
  if fid==-1
    error(['File  ',outputdir,'  not found.']);
  end
  consv = fscanf(fid,'%e',[mcons+1 inf]);
  status = fclose(fid);

  t  = consv(1,:)';
  nt = length(t);

  % line for plots %
  gamma = -0.1533;
  tsmall = t(1: min(find( t > 15 ) ) );
  y = consv(2,1) * exp(gamma*tsmall);

  % Electric Field
  normE = consv(2,:)';

  figure(1);

  q = semilogy(t,normE,'b-', 'linewidth', 3);
  hold off;

  % L1 norm of f
  figure(2);

  axis on; box on; grid off;
  axis([0 60 -1e-12 1e-12]);
  qc   = consv(3,:)';
  dqc  = qc - qc(1);
  hold on;
  p = plot(t,dqc,'b-', 'linewidth', 3);
  hold on;
  q = plot(t, zeros(size(t)), 'k-' );
  set(q, 'linewidth', 3 );
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
  p = plot(t,dqc,'b-');
  set(p, 'linewidth', 3 );
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
  axis([0 60 -0.01 1]);
  hold off;

  for n=1:5
    figure(n);
    set(gca,'fontsize',16);
    l = legend('SL2', 'SL4');
    set(l, 'location', 'northwest');
  end
  figure(1);
  l = legend('SL2', 'SL4');
  set(l, 'location', 'north');
  figure(3);
  l = legend('SL2', 'SL4');
  set(l, 'location', 'southwest');

  % print all the pretty pictures
  print(1, '-depsc2', 'Strong_128_ElectricField.eps' );
  print(2, '-depsc2', 'Strong_128_L1f.eps' );
  print(3, '-depsc2', 'Strong_128_L2f.eps' );
  print(4, '-depsc2', 'Strong_128_Energy.eps' );
  print(5, '-depsc2', 'Strong_128_Entropy.eps' );

end
