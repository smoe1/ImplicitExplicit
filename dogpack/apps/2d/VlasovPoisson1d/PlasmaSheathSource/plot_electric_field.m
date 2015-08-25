function plot_electric_field(points_per_dir,outputdir_in, point_type )
%PLOT_ELECTRIC_FIELD
%
%  This is essentially a copy of plotdog1, but used to plot the electric field
%  instead.
%
%  Usage: plotdog1( points_per_dir, outputdir_in, point_type  )
%
% points_per_dir = number of points used per cell for plotting.  Default = 1
%                  witout any arguments supplied 
%
% outputdir_in = location for output directory.  defulat = 'output'
%
% point_type = 1 : Points are forced to be uniformly spaced.
% point_type = 2 : Gaussian Quadrature points are used
%

  global outputdir
  outputdir='output';
  
  if(nargin<1)
    points_per_dir = 1;
  end 

  if (ischar(points_per_dir))
    points_per_dir = str2num(points_per_dir);
  end

  if(nargin>1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end

  if points_per_dir<1
    points_per_dir = 1;
  end
  
  if (nargin>2)
    if(ischar(point_type))
      point_type = str2num(point_type);
    end
  else
    point_type=1;
  end
  
  if point_type~=1 && point_type~=2
    point_type = 1;
  end  
  
  if point_type==2 && points_per_dir>5
    disp('too many points chosen, using 5 points instead');
    points_per_dir=5;
  end


  disp(' ');
  disp(['       outputdir = ',outputdir]);
  disp(['  points_per_dir = ',num2str(points_per_dir)]);
  %disp(' ');

  fids  = fopen([outputdir,'/qhelp1d.dat'],'r');
  if fids==-1
  error(['File  ',outputdir,'/qhelp1d.dat  not found.']);
  end
  nplot = fscanf(fids,'%d',1);  % number of plots
  meqn  = fscanf(fids,'%d',1);  % number of equations
  maux  = fscanf(fids,'%d',1);  % number of auxiliary equations >= 0
  meth1 = fscanf(fids,'%d',1);  % order of accuracy 
  mx    = fscanf(fids,'%d',1);  % number of grid elements
  xlow  = fscanf(fids,'%e',1);  % left hand side of domain
  xhigh = fscanf(fids,'%e',1);  % right hand side of domain
  dx    = fscanf(fids,'%e',1);  % grid width (this should be consistent with above ...)
  fclose(fids);

  disp(['          melems = ', int2str(mx)]);
  disp(' ');

  % Grid information
  mx_old = mx;
  mx = mx*points_per_dir;
  dx_old = dx;
  dx = (xhigh-xlow)/mx;

  % TODO - this needs to be modified ...
  if( point_type == 1 )

      xc = transpose(linspace(xlow+dx/2,xhigh-dx/2,mx));
      xc_old = transpose(linspace(xlow+dx_old/2,xhigh-dx_old/2,mx_old));

      % Linearly spaced points
      s1d    = -1.0 + (2.0*(1:points_per_dir)-1.0)/points_per_dir;

  else

    sq3 = sqrt(3);
    sq5 = sqrt(5);
    sq7 = sqrt(7);
   
    % quadrautre points (TODO - add in 6th order story ... )
    if (points_per_dir==1)
      s1d = 0.0;
    elseif (points_per_dir==2)
      s1d = [-1.0/sq3; 1.0/sq3];
    elseif (points_per_dir==3)
      s1d = [-sq3/sq5; 0.0e0; sq3/sq5];
    elseif (points_per_dir==4)
      s1d = [-sqrt(3.0+sqrt(4.8))/sq7; -sqrt(3.0-sqrt(4.8))/sq7; ...
              sqrt(3.0-sqrt(4.8))/sq7;  sqrt(3.0+sqrt(4.8))/sq7];
    elseif (points_per_dir==5)
      s1d = [-sqrt(5.0 + sqrt(40.0/7.0))/3.0; ...
             -sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
              0.0; ...
              sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
              sqrt(5.0 + sqrt(40.0/7.0))/3.0];
    end

    xc = zeros(mx_old*points_per_dir,1);
        
    xline = (xlow+dx_old/2):dx_old:(xhigh-dx_old/2);

    kk=1;
    for i=1:mx_old                 
      xc( kk:(kk+points_per_dir-1) ) = xline(i)+(dx_old/2)*s1d(:);
      kk=kk+points_per_dir;
    end

  end


  % Sample basis functions on mesh
  phi = GetCart1Legendre(meth1, s1d);

  q=-1;

  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
    disp(' ')
if isempty(m)
    m=1;
    end

    kn = 0;

    n = 0;
    nf = 0;
    n1 = -1;

  while(nf~=-1)
        nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
                ' ) [type -1 or q to quit] ? ']);
    if isempty(nf)
      n1 = n1 + 1;
      nf = 0;
      else
      n1 = nf;
      end
      if n1> nplot
      disp(' ');
      disp(' End of plots ');
      disp(' ');
      n1 = nplot;
      end
  
    if (nf~=-1)
  
      %% Solution -- q
      % solution should be found in file
      %     outputdir/q[n1].dat
      fname = [outputdir,'/',num2str(n1+10000),'.dat'];
  
      % replace the 1000's digit by the letter q
      fname(length(outputdir)+2) = 'e';
  
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
  
      %% Aux variables -- aux
      if (maux>0)
          fname(length(outputdir)+2) = 'a';
          fids = fopen(fname,'r');
          if fids==-1
              error(['File  ',fname,'  not found.']);
          end
          time = fscanf(fids,'%e',1);
          atmp = fscanf(fids,'%e',[1,inf]);
          fclose(fids);
          atmp = transpose(atmp);
          acoeffs  = reshape(atmp,mx_old,maux,meth1);
          clear atmp;
          aux = zeros(mx,maux);
          for i=1:mx_old
              for ma=1:maux
                  for ii=1:points_per_dir              
                      v1(1:meth1,1) = phi(ii,:);
                      v2(1:meth1,1) = acoeffs(i,ma,:);
                      aux((i-1)*points_per_dir+ii,ma) = transpose(v1)*v2;
                  end
              end
          end
      end
  
      % USER SUPPLIED FUNCTION
      plotq1;        
      end
  
   end
   disp(' ')
