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
 
  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
    ' ) ? ']);
  disp(' ')
  if isempty(m)
    m=1;
  end

  % INFO --------------------------------------------------------------

  mcons = 5;
  fid = fopen([outputdir '/Efield.dat']);
  consv = fscanf(fid,'%e',[mcons+1 inf]);
  status = fclose(fid);

t  = consv(1,:)';
nt = length(t);

% line for plots %
gamma = -0.1533;
y = consv(2,1) * exp(gamma*t);

% Electric Field
normE = consv(2,:)';

figure(1);

clf
set(gca,'fontsize',16);
t1=title([ 'Weak Landau Damping: Electric Field vs. Time']);
set(t1,'fontsize',16);
axis on; box on; grid off;
hold off

q = semilogy(t,normE,'b-');
set(q, 'LineWidth', 3 );
hold on;
p = semilogy(t,y,'k-');
set(p, 'LineWidth', 3 );

figure(2);
clf

% L1 norm of f
subplot(2,2,1);

%set(gca,'fontsize',16);
%t1=title([ 'Weak Landau Damping: Conserved Quantities vs. Time']);
%set(t1,'fontsize',16);
%set(gca,'fontsize',16);

axis on; box on; grid off;
%axis([0 60 -0.001 0.001]);
qc   = consv(3,:)';
dqc  = qc - qc(1);
p = plot(t,dqc);
set(p, 'linewidth', 3 );
hold on;
q = plot(t, zeros(size(t)), 'k-' );
set(q, 'linewidth', 3 );
hold off;


% L2 norm of f
subplot(2,2,2);
qc   = consv(4,:)';
dqc  = qc - qc(1);
p = plot(t,dqc);
set(p, 'linewidth', 3 );
hold on;
q = plot(t, zeros(size(t)), 'k-' );
set(q, 'linewidth', 3 );
hold off;

% Total Energy
subplot(2,2,3);
qc   = consv(5,:)';
dqc  = qc - qc(1);
p = plot(t,dqc);
set(p, 'linewidth', 3 );
hold on;
q = plot(t, zeros(size(t)), 'k-' );
set(q, 'linewidth', 3 );
hold off;

% Entropy 
subplot(2,2,4);
qc   = consv(6,:)';
dqc  = qc - qc(1);
p = plot(t, -dqc);
set(p, 'linewidth', 3 );
hold on;
q = plot(t, zeros(size(t)), 'k-' );
set(q, 'linewidth', 3 );
hold off;

figure(1);

print(1, '-depsc2', 'WeakHSL4_64_128_ElectricField.eps' );
