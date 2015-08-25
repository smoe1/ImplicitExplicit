function [time,E,E2,aux,dx] = plotE2(outputdir_in)

  global outputdir
  outputdir='output';
  
  if(nargin>0)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end

  format long e;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  FIND OUT IF CARTESIAN OR UNSTRUCTURED GRID
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fids  = fopen([outputdir,'/qhelp.dat'],'r');
  if fids==-1
    error(['File  ',outputdir,'/qhelp.dat  not found.']);
  end
  GridType = fread(fids, 9, 'char')';
  if (GridType=='Cartesian')
    GridType='Cartesian   ';
  elseif (GridType=='Unstructu')
    GridType='Unstructured';
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(' ');
  disp(['        GridType = ',GridType]);
  disp(['       outputdir = ',outputdir]);
  disp(' ');

  if (GridType=='Cartesian   ')
    
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
  
    % Grid information
    kmax = get_kmax(meth1);
    dx = (xhigh-xlow)/mx;
    dy = (yhigh-ylow)/my;
    
    E    = zeros(nplot,mx,1);
    E2   = zeros(nplot,1);
    time = zeros(nplot,1);

    for n1 = 1:nplot
     
      %% Aux variables -- aux

      %%% cutting out guts from read_state2_cart %%%
      num_frame = n1;
      varname = 'a';
      basefilename = [outputdir '/' varname sprintf('%04d',num_frame)];
      filename = [basefilename '.dat'];
      fids = fopen(filename,'r');
      assert(fids~=-1, ['could not open file ' filename]);
      time(n1) = fscanf(fids,'%e',1);
      aux_data = fscanf(fids,'%e',[inf]);
      fclose(fids);

      dims = [mx,my,maux,kmax];
      aux = reshape(aux_data,dims);
%     permutation=[1,2,4,3];
%     aux = permute(aux,permutation);

      E(n1,:,1) = aux(:,1,2,1);
      E2(n1) = sqrt( sum( aux(:,9,2,1) .* aux(:,9,2,1) ));

    end

  end
  disp(' ')

  figure(1);
  clf
  plot(time,E2);
  set(gca,'fontsize',16);
  t1 = title(['Landau Damping \alpha = 0.5,  mx = ', int2str(mx), '     [DoGPack]']);
  set(t1,'fontsize',16);
  xlabel('time','FontSize',16) 
  ylabel('|E|_2','FontSize',16)

end
