function [xc,qsoln] = get_mhd(n1,points_per_dir,outputdir)
    
    fids  = fopen([outputdir,'/qhelp.dat'],'r');
    if fids==-1
        error(['File  ',outputdir,'/qhelp.dat  not found.']);
    end
    nplot = fscanf(fids,'%d',1);
    meqn  = fscanf(fids,'%d',1);
    maux  = fscanf(fids,'%d',1);
    meth1 = fscanf(fids,'%d',1);
    mx    = fscanf(fids,'%d',1);
    xlow  = fscanf(fids,'%e',1);
    xhigh = fscanf(fids,'%e',1);
    dx    = fscanf(fids,'%e',1);
    fclose(fids);
    
    % Grid information
    mx_old = mx;
    mx = mx*points_per_dir;
    dx_old = dx;
    dx = (xhigh-xlow)/mx;
    xc = transpose(linspace(xlow+dx/2,xhigh-dx/2,mx));
    xc_old = transpose(linspace(xlow+dx_old/2,xhigh-dx_old/2,mx_old));
    
    % Sample basis functions on mesh
    % size of phi = (points_per_dir,meth1);
    phi = SampleBasis1(points_per_dir,meth1);
    
    %% Solution -- q
    % solution should be found in file
    %     outputdir/q[n1].dat
    fname = [outputdir,'/',num2str(n1+10000),'.dat'];
    
    % replace the 1000's digit by the letter q
    fname(length(outputdir)+2) = 'q';
    
    fids = fopen(fname,'r');
    if fids==-1
        error(['File  ',fname,'  not found.']);
    end
    
    time = fscanf(fids,'%e',1);
    qtmp = fscanf(fids,'%e',[1,inf]);
    fclose(fids);
    qtmp = transpose(qtmp);
    qcoeffs  = reshape(qtmp,mx_old,meqn,meth1);
    clear qtmp;
    qsoln = zeros(mx,meqn);
    for i=1:mx_old
        for me=1:meqn
            for ii=1:points_per_dir
                v1(1:meth1,1) = phi(ii,:);
                v2(1:meth1,1) = qcoeffs(i,me,:);
                qsoln((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
            end
        end
    end