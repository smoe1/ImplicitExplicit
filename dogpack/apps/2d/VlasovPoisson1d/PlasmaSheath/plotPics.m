function plot_E( mxsize )

  for n=1:5
    figure(n);
    clf;
  end

  mcons = 5;  % number of conserved quantities 

% fid = fopen(['outputSL4_' int2str(mxsize) '/Efield.dat']);

% if fid==-1
%   error(['File  ',outputdir,'  not found.']);
% end
% consv = fscanf(fid,'%e',[mcons+1 inf]);
% status = fclose(fid);

% t  = consv(1,:)';
% nt = length(t);

% % Electric Field
% normE = consv(2,:)';

% figure(1);

% q = semilogy(t,normE,'ro');
% hold on;

% % L1 norm of f
% figure(2);

% qc   = consv(3,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'r--');
% set(p, 'linewidth', 3 );
% hold off;

% % L2 norm of f
% figure(3);
% qc   = consv(4,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'r--');
% set(p, 'linewidth', 3 );
% hold off;

% % Total Energy
% figure(4);
% qc   = consv(5,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'r--');
% set(p, 'linewidth', 3 );
% hold off;

% % Entropy 
% figure(5);
% qc   = consv(6,:)';
% dqc  = qc - qc(1);
% p = plot(t, dqc,'r--');
% set(p, 'linewidth', 3 );
% hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Data Set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid = fopen([strcat('outputSL2_', int2str(mxsize)) '/Efield.dat']);
% if fid==-1
%   error(['outputdir not found.']);
% end
% consv = fscanf(fid,'%e',[mcons+1 inf]);
% status = fclose(fid);

% t  = consv(1,:)';
% nt = length(t);

% % Electric Field
% normE = consv(2,:)';

% figure(1);

% q = semilogy(t,normE,'go');
% hold on;

% % L1 norm of f
% figure(2);

% qc   = consv(3,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'g--');
% set(p, 'linewidth', 3 );
% hold off;

% % L2 norm of f
% figure(3);
% qc   = consv(4,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'g--');
% set(p, 'linewidth', 3 );
% hold off;

% % Total Energy
% figure(4);
% qc   = consv(5,:)';
% dqc  = qc - qc(1);
% p = plot(t,dqc,'g--');
% set(p, 'linewidth', 3 );
% hold off;

% % Entropy 
% figure(5);
% qc   = consv(6,:)';
% dqc  = qc - qc(1);
% p = plot(t, dqc,'g--');
% set(p, 'linewidth', 3 );
% hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%% Third data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

  % plot the first line:
% if mxsize == 128
%   x1 = [2.4470588235294137e+00   1.1140499778281017e+00];
%   x2 = [4.5490196078431318e+00   6.0299329769456655e-01];
% elseif mxsize == 64
%   x1 = [2.4566929133858277e+00   1.1137239792298530e+00];
%   x2 = [4.5354330708661390e+00   6.0298185286541339e-01];
% else
%   disp(['  warning, no data for plotting correct lines given']);
%   x1 = [0 1];
%   x2 = [1 1];
% end

% gamma = (log(x2(2)) - log(x1(2))) / (x2(1) - x1(1) )
% t0 = -log( x1(2) ) / gamma + x1(1);

% tsmall = 0:0.1:9;
% y  = exp( gamma*(tsmall-t0) );
% p = semilogy(tsmall,y,'k-', 'linewidth', 3);

  % plot the second line
% if mxsize == 128
%   x1 = [2.3278431372548948e+01   1.1330961736590262e-01];
%   x2 = [2.8423529411764616e+01   1.7231616364197944e-01];
% elseif mxsize == 64
%   x1 = [2.3307086614173194e+01   1.1387731193552551e-01];
%   x2 = [2.8409448818897594e+01   1.7514508832396741e-01];
% else
%   disp(['  warning, no data for plotting correct lines given']);
%   x1 = [0 1];
%   x2 = [1 1];
% end

% gamma2 = (log(x2(2)) - log(x1(2))) / (x2(1) - x1(1) )
% t0 = -log( x1(2) ) / gamma2 + x1(1);

% tsmall = 21:0.1:38;
% y2 = exp( gamma2*(tsmall-t0) );
% hold on;
% plot(tsmall, y2, 'k-', 'linewidth', 3);

% hold off;

  % L1 norm of f
  figure(2);

  axis([0 60 -1e-12 1e-12]);
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
% t1=title([ 'L_2(Electric Field)  \gamma_1 = ', ...
%   sprintf('%2.4f', gamma), ',  \gamma_2 = ', sprintf('%2.4f', gamma2)] );
  set(t1,'fontsize',16);
  axis([0 60 1e-4 1e2]);
  l = legend('SL2', 'SL4', 'location', 'north');

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
  axis([0 60 -0.1 0.02]);
  l = legend('SL2', 'SL4', 'location', 'southwest');

  % energy
  figure(4);
  set(gca,'fontsize',16);
  t1=title(['Energy( t ) - Energy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 60 -0.01 0.04]);
  l = legend('SL2', 'SL4', 'location', 'northwest');

  % entropy
  figure(5);
  set(gca,'fontsize',16);
  t1=title(['Entropy( t ) - Entropy( t=0 )']);
  set(t1,'fontsize',16);
  axis([0 60 -0.1 1]);
  l = legend('SL2', 'SL4', 'location', 'northwest');
  set(l, 'fontsize', 16 );

  % print all the pretty pictures
% name1 = strcat( 'Strong_', strcat( int2str(mxsize), '_ElectricField.eps' ) );
% name2 = strcat( 'Strong_', strcat( int2str(mxsize), '_L1f.eps' ) );
% name3 = strcat( 'Strong_', strcat( int2str(mxsize), '_L2f.eps' ) );
% name4 = strcat( 'Strong_', strcat( int2str(mxsize), '_Energy.eps' ) );
% name5 = strcat( 'Strong_', strcat( int2str(mxsize), '_Entropy.eps' ) );
% print(1, '-depsc2', name1 );
% print(2, '-depsc2', name2 );
% print(3, '-depsc2', name3 );
% print(4, '-depsc2', name4 );
% print(5, '-depsc2', name5 );

end
