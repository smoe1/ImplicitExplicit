%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
clf;
pz=plot(xc,qsoln(:,m),'r-');
set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
%axis([0 1 -0.6 1.6]);
axis([0 1 -1.1+0.5 1.1+0.5]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);


%sp = zeros(mx_old,1);
%se = zeros(mx_old,1);
%flag = zeros(mx_old,1);
%for i=1:mx_old
%  Q2sum = 0.0;
%      
%  for k=1:(meth1-2)
%    Q2sum = Q2sum + qcoeffs(i,1,k)^2;
%  end
%      
%  Q2kmax = qcoeffs(i,1,meth1-1)^2;
%  Q2sum = Q2sum + Q2kmax;
%
%  if Q2kmax>1.0e-10
%    sp(i) = (Q2kmax/Q2sum);
%  else
%    sp(i) = 0.0;
%  end
%  
%  Q2kmax = qcoeffs(i,1,meth1)^2;
%  Q2sum = Q2sum + Q2kmax;
%  
%  if Q2kmax>1.0e-10
%    se(i) = (Q2kmax/Q2sum);
%  else
%    se(i) = 0.0;
%  end
%  
%  if (se(i)>=(0.1*(meth1)^(-4)))
%    flag(i) = 1;
%  end
%  
%end

%figure(2);
%clf;
%pz=plot(xc_old,se,'bo');
%%hold on
%%pzz=plot(xc_old,sp,'rx');
%hold off;
%axis on; box on; grid off;
%set(gca,'plotboxaspectratio',[2 1 1]);
%set(gca,'fontsize',16);

%figure(3);
%clf;
%pz=plot(xc_old,flag,'bo');
%hold off;
%axis on; box on; grid off;
%set(gca,'plotboxaspectratio',[2 1 1]);
%set(gca,'fontsize',16);

%% solve for the exact solution: %%
%t = 0.15;
%t = time;
%q0  = @(z)(sin(2*pi*z));
%qp = @(z)(2*pi*cos(2*pi*z));
%f  = @(z)(q0(z)*t + z - xc);
%fp = @(z)(qp(z)*t + 1);

%if( abs(t-0.10) < 1e-10 )
%  disp('  computing exact solution');
%  xi = xc;   % initial guess
%  tol = 1e-15;
%  err = max(abs( f(xi) ) );
%  num_iter = 0;
%  while( err > tol )
%  	xi = xi - (f(xi) ./ fp(xi));
%  	err = max(abs(f(xi)));
%	num_iter = num_iter+1;
%	if( num_iter > 100 )
%		disp('   too many iterations' );
%		break;
%	end
%  end
%
%  qex = q0(xc-t*q0(xi));
%  disp(['   max(qex-q(xi)) = ', num2str(max(qex-q(xi)),'%0.8e')]);
%
%  hold on;
%  pt = plot(xc,qex, 'r-');
%  set(pt,'linewidth',1);

%  disp(['   min(xi) = ', num2str(min(xi)), '    max(xi) = ', num2str(max(xi))]);

  %  print error in exact vs computed solution
%  err2 = norm(qsoln-qex,2)/norm(qex,2);
%  disp(' ');
%  disp(['   dx = ',num2str(dx,'%0.8e'),['         err2 = ' ...
%                      ''],num2str(err2,'%0.8e')]);
%  disp(' ');

%end
