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
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['Height at t = ',num2str(time),'     [DoGPack]']);
axis([0 1 0.75 3.25]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:1:3);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

figure(2);
clf;  
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['hu^1(x,t) at t = ',num2str(time),'     [DoGPack]']);  
axis([0 1 -1.6 1.6]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1.5:0.5:1.5);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);


if( abs(time-0.1) < 1e-5 )
    % read in data from exact solution (highly resolved run) and compare with
	% current results.  should run a highly resolved computation and output
	% results to folder output_ex in order for this to work
	%fids  = fopen([outputdir,'/qhelp.dat'],'r');
	fids  = fopen(['output_ex/qhelp.dat'],'r');
	if fids==-1
	  error(['File  ','output_ex/qhelp.dat  not found.']);
	end
	nplot_ex = fscanf(fids,'%d',1);
	meqn_ex  = fscanf(fids,'%d',1);
	maux_ex  = fscanf(fids,'%d',1);
	meth1_ex = fscanf(fids,'%d',1);
	mx_ex    = fscanf(fids,'%d',1);
	xlow_ex  = fscanf(fids,'%e',1);
	xhigh_ex = fscanf(fids,'%e',1);
	dx_ex    = fscanf(fids,'%e',1);
	fclose(fids);
	xc_ex = transpose(linspace(xlow+dx/2,xhigh-dx/2,mx_ex));

    %% Exact Solution -- q
%    fname = ['output_ex/',num2str(n1+10000),'.dat'];
%	size(fname)
%    fname(8) = 'q_ex';
%    fids = fopen(fname,'r');

% hack:  Isample should be opening the 'last' file printed.  THIS WILL ONLY WORK
% WHEN USER USES EXACTLY 10 output files!
    fids = fopen('output_ex/q0010.dat');
    time = fscanf(fids,'%e',1);
    qtmp_ex = fscanf(fids,'%e',[1,inf]);
    fclose(fids);
    qtmp_ex = transpose(qtmp_ex);
%    qcoeffs_ex  = reshape(qtmp_ex,mx_old_ex,meqn_ex,meth1);
    qcoeffs_ex  = reshape(qtmp_ex,mx_ex,meqn,meth1_ex);
    clear qtmp;
    qsoln_ex = zeros(mx_ex,meqn_ex);


    phi_ex = SampleBasis1(points_per_dir,meth1);
    %for i=1:mx_old
    for i=1:mx_ex
      for me=1:meqn
        for ii=1:points_per_dir
          v1(1:meth1_ex,1) = phi_ex(ii,:);
          v2(1:meth1_ex,1) = qcoeffs_ex(i,me,:);
          qsoln_ex((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
        end
      end
    end

   ratio = length(qsoln_ex) / length(qsoln);
   if( mod(ratio,3) ~= 0)
       disp('   bad ratio between qex and qsoln');
	   disp(['     ratio = ',  int2str(ratio) ] );
	   return;
   end

   % index of first grid cell the two solutions have in common
   nstart = int32 ((ratio + 1) / 2);

   % sample grid points
   Isample = int32( ratio * (0:(mx-1)) ) + int32(nstart);

%   disp(['            size(Isample) = ', int2str(size(Isample)) ]);
%   disp(['        size(qsoln) = ', int2str(size(qsoln))]);
%   disp(['     size(qsoln_ex) = ', int2str(size(qsoln_ex))]);
%   disp([' size(qsoln_ex(Isample,:) = ', int2str(size(qsoln_ex(Isample,:)))]);

   err = norm( (qsoln_ex(Isample,1) - qsoln(:,1)),2) ./ norm(qsoln_ex(Isample,1), 2);
   disp(['       err = ', num2str(err, '%0.8e')]);

   % plot the 'exact solution' on top of this solution
   figure(1)
   hold on
   pz=plot(xc,qsoln_ex(Isample,1),'-r');

   figure(2)
   hold on
   pz=plot(xc,qsoln_ex(Isample,2),'-r');

end % end of comparing solution with highly resolved run
