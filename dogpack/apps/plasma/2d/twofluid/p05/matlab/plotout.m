function plotout(outputdir_in)

  setarg_outputdir;
  
  all_names = { ...
     ... % state variables
     'rho_i', 'Mi', 'Ni', ...
     'B', 'E', 'Psi', ...
     'bogus', ...
      ... % quantities to be constructed
     'sgm', 'dvB', 'dvE', 'u', 'J', 'Ji', ...
     'p_i', 'p', ...
     ... % Ohm's law terms
     'Bxu', 'Hal', 'pek', 'J_t', 'DFJ', 'dtJ', 'Eck', ...
     ... % Ohm's law--related terms
     'Mi_t', ...
  };

  plotout2(outputdir, all_names);
end

function plotout2(outputdir, all_names)

  format long e;

  % TITLE INFO --------------------------------------------------------------
  disp(' ')
  disp(['DoGPack: plotting output from ' outputdir]);
  disp(' ')
  % END TITLE INFO ----------------------------------------------------------
  
  % get configuration parameters and mesh data
  %
  params = get_parameters(outputdir);
  %
  % aliases for fields
  %
  space_order = params.space_order;
  mx        = params.mx;
  my        = params.my;
  xlow      = params.xlow;
  xhigh     = params.xhigh;
  ylow      = params.ylow;
  yhigh     = params.yhigh;
  datafmt   = params.datafmt;
  %
  BCs                = params.BCs;
  domain_scaling     = params.domain_scaling;
  enforced_symmetry  = params.enforced_symmetry;
  mass_ratio         = params.mass_ratio;
  cc                 = params.cc;
  cs_light           = params.cs_light;
  %
  one_over_epsilon = params.one_over_epsilon;
  ion_mass = params.ion_mass;
  elc_mass = params.elc_mass;

  [xl,yl,xc,yc]=get_grid2(params);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  allNames = ComponentNames(all_names);
  index = get_name2idx(allNames);
  %
  % define global variable component names
  % (deprecated)
  for i=1:length(all_names)
    %eval([char(all_names(i)) '=' '''' char(all_names(i)) ''';']);
    eval(['global ' char(all_names(i)) ';'])
    eval([char(all_names(i)) '=' num2str(i) ';'])
  end
  %
  all_num_components=zeros(1,length(all_names));
  all_num_components(index.rho_i   ) = 1;
  all_num_components(index.Mi      ) = 3;
  all_num_components(index.Ni      ) = 1;
  all_num_components(index.B       ) = 2;
  all_num_components(index.E       ) = 1;
  all_num_components(index.Psi     ) = 1;
  %
  all_num_components(index.Bxu) = 1;
  all_num_components(index.Hal) = 1;
  all_num_components(index.pek) = 1; % change to 2?
  all_num_components(index.J_t) = 1;
  all_num_components(index.DFJ) = 1;
  all_num_components(index.dtJ) = 1;
  all_num_components(index.Eck) = 1;
  all_num_components(index.sgm) = 1;
  all_num_components(index.dvB) = 1;
  all_num_components(index.dvE) = 1;
  all_num_components(index.u  ) = 3;
  all_num_components(index.J  ) = 3;
  all_num_components(index.Ji ) = 3;
  all_num_components(index.p_i) = 1;
  all_num_components(index.p  ) = 1;
  all_num_components(index.Mi_t)= 3;
  %
  all_flip=ones(1,length(all_names));
  all_flip([index.B,index.Psi,index.dvB]) = -1; % flip pseudovector quantities
  %
  state_start_components=zeros(1,index.bogus);
  state_start_components(1) = 1;
  for i=2:index.bogus
    state_start_components(i) = ...
       state_start_components(i-1)+all_num_components(i-1);
  end
  last_state_idx = index.bogus-1;

  % deprecated
  %
  global state_components;
  state_components=cell(1,last_state_idx);
  for i=1:last_state_idx
    state_components{i}=state_start_components(i):state_start_components(i+1)-1;
  end
  assert(params.meqn==state_components{last_state_idx}(end));

  state_idx = struct;
  for i=1:last_state_idx
    state_idx.(all_names{i}) = state_start_components(i):state_start_components(i+1)-1;
  end

  %state_components
  %
  % 0 = 4 = do not display vector field
  % 1 = represent vector field with contours
  % 2 = represent vector field with vectors
  % 3 = represent vector field both ways
  all_noDiv=ones(1,length(all_names))*2;
  %for i=1:length(all_names)
  %  if(all_num_components(i)==3)
  %    all_noDiv(i) = 2;
  %  end
  %end
  all_noDiv(B) = 1;
  %
  all_color_ranges=zeros(2,length(all_names));
  all_color_ranges(:,index.rho_i   ) = [0,1];
  all_color_ranges(:,index.Mi      ) = [-0,0]; %[-0.2,.2];
  all_color_ranges(:,index.Ni      ) = [0,1];
  all_color_ranges(:,index.B       ) = [-.2,.2];
  all_color_ranges(:,index.E       ) = [-.2,.2];
  all_color_ranges(:,index.Psi     ) = [-.05,.05];
  all_color_ranges(:,index.Bxu     ) = [-.2,.2];
  all_color_ranges(:,index.Hal     ) = [-.2,.2];
  all_color_ranges(:,index.pek     ) = [-.2,.2];
  all_color_ranges(:,index.J_t     ) = [-.2,.2];
  all_color_ranges(:,index.DFJ     ) = [-.2,.2];
  all_color_ranges(:,index.dtJ     ) = [-.2,.2];
  all_color_ranges(:,index.Eck     ) = [-.2,.2];
  all_color_ranges(:,index.sgm     ) = [-1,1]*2e-1;
  all_color_ranges(:,index.dvB     ) = [-.03,.03];
  all_color_ranges(:,index.dvE     ) = [-10,10];  % why is this so large?
  all_color_ranges(:,index.u       ) = [-0,0];
  all_color_ranges(:,index.J       ) = [-0,0];
  all_color_ranges(:,index.Ji      ) = [-0,0];
  all_color_ranges(:,index.p       ) = [ 0,1];

  all_vector_scales=zeros(1,length(all_names));
  all_vector_scales(:,index.rho_i   ) = 0;
  all_vector_scales(:,index.Mi      ) = all_color_ranges(2,index.Mi);
  all_vector_scales(:,index.Ni      ) = 0;
  all_vector_scales(:,index.B       ) = 1;
  all_vector_scales(:,index.E       ) = all_color_ranges(2,index.E);
  all_vector_scales(:,index.Psi     ) = 0;
  all_vector_scales(:,index.Bxu     ) = all_color_ranges(2,index.Bxu);
  all_vector_scales(:,index.Hal     ) = all_color_ranges(2,index.Hal);
  all_vector_scales(:,index.pek     ) = all_color_ranges(2,index.pek);
  all_vector_scales(:,index.J_t     ) = all_color_ranges(2,index.J_t);
  all_vector_scales(:,index.DFJ     ) = all_color_ranges(2,index.DFJ);
  all_vector_scales(:,index.dtJ     ) = all_color_ranges(2,index.dtJ);
  all_vector_scales(:,index.Eck     ) = all_color_ranges(2,index.Eck);
  all_vector_scales(:,index.sgm     ) = 0;
  all_vector_scales(:,index.dvB     ) = 0;
  all_vector_scales(:,index.dvE     ) = 0;
  all_vector_scales(:,index.u       ) = all_color_ranges(2,index.u);
  all_vector_scales(:,index.J       ) = all_color_ranges(2,index.J);
  % so that user can enter q to quit
  q=-1;
  qq=-2;

  mf_display=0;
  while(mf_display~=-1)
    % display choices
    for i=1:length(all_names)
      disp([num2str(i), ' = ', char(all_names(i))]);
    end
    
    default_input_string='[E,Eck,Bxu;Hal,dtJ,pek]''';
    input_string = input([ 'Which component[s] of q do you want to plot ' ...
            default_input_string ' ) ? '],'s');
    disp(' ')
    if(isempty(input_string)||isempty(mf_display))
      % default components
      %mf_display=[E,Eck,Bxu;Hal,dtJ,pek]';
      mf_display = eval(default_input_string)
    else
      mf_display = eval(input_string);
    end
    %
    if(length(mf_display)==1 && (mf_display==q || mf_display==qq))
      break;
    end
    %
    % identify entries with '0' as bogus, so a blank
    % spot can be put there.
    indexes_of_invalid_mf_entries = ...
       find(mf_display<1 | mf_display> length(all_names));
    mf_display(indexes_of_invalid_mf_entries)=index.bogus;
    %mf_display_files=mf_display
    %indexes_of_nonfile_mf_entries = find(mf_display>length(all_fname_prefixes));
    %mf_display_files(indexes_of_nonfile_mf_entries)=bogus;
    %%mf_display_isfile=(mf_display <= length(all_fname_prefixes));

    % restrict to selected subset of arrays
    %
    num_components = all_num_components (mf_display);
    names          = all_names          (mf_display);
    flips          = all_flip           (mf_display);
    noDiv          = all_noDiv          (mf_display);
    color_ranges   = all_color_ranges   (:,mf_display);
    vector_scales  = all_vector_scales  (:,mf_display);
    %fname_prefixes = all_fname_prefixes (mf_display_files);

    % create/raise the plot windows to be used
    use_subplots=0;
    pane_num=0;
    if(size(mf_display,1)>1)
      use_subplots=1;
    end
    if(use_subplots)
      pane_num = length(all_names)+1;
      figure(pane_num);
    else % make array of subplots
      for i=1:numel(mf_display)
        figure(mf_display(i));
        % leave the contents of the figure untouched
        % in case the user wants to use this feature
        % to make a comparison.
        %clf;
      end
    end

    nf = 0;
    n1 = -1;
    while(nf~=-1)
      nf  = input([ ' Plot which frame ( 0 - ',num2str(params.nout),...
            ' ) [type -1 or q to quit] ? ']);
      if (nf==q)
        break
      end
      if(nf==qq)
        return
      end
      if isempty(nf)
        nf = n1 + 1;
      end
      if max(nf)> params.nout
        disp(' ');
        disp(' End of plots ');
        disp(' ');
        nf = params.nout;
      end

      % in case the user entered something like 0:5:40 (every fifth frame)
      next_frame=0;
      for idx=1:length(nf)
        n1=nf(idx);
        if(next_frame==0)
          next_frame=1;
        else
          go_next_frame  = input([  ' Plot frame ' num2str(n1) ' (or q)?']);
          if (go_next_frame==qq)
            return;
          elseif (go_next_frame==q)
            break
          end
        end
        %
        READ_FROM_RESTART=1;
        % in matlab it seems faster just to read the whole output
        % than to work with low-level access to array slices;
        % in any case it seems that most of the time is
        % consumed by plotting rather than reading or reconstructing.
        % So the procedure is:
        % 1. read everything,
        % 2. sample and plot the relevant components
        READ_ALL_AT_ONCE=1;
        basefilename = [outputdir sprintf('/q%04d',n1)];
        [state,time]=read_state2_cart(datafmt, outputdir, n1, 'q',...
          mx, my, params.meqn, get_kmax(space_order), 1:params.meqn);
        basename = sprintf('%04d',n1);
        for i=1:numel(mf_display)
           if(mf_display(i)==bogus) % bogus indicates an empty subplot
              %subplot(size(mf_display,1),size(mf_display,2),i);
              %clf; % need something that will only clear the subplot
           else

              disp(['calculating term ' names{i}]);
              tic;
              % see compute_dtJ for documentation on data representation
              % a better way to handle this would be:
              % + create a "class" (cell array) which contains
              %   - a data format indicator field and
              %   - the array
              % + convert format as needed
              % + pass to each compute_* function only the data it requires
              % + unix philosophy: split up into small routines each of which
              %   does one thing well.
              %
              if(mf_display(i)<=last_state_idx)
                components=state_components{mf_display(i)};
                data=sample_state2(state(:,:,:,components),params.space_order);
              else
                if(mf_display(i)==sgm)
                  data = compute_sgm(state);
                elseif(mf_display(i)==u)
                  data = compute_u(state);
                elseif(mf_display(i)==J)
                  data = compute_J(state);
                elseif(mf_display(i)==Ji)
                  data = compute_Ji(state);
                elseif(mf_display(i)==dvB)
                  data = compute_dvB(state);
                elseif(mf_display(i)==dvE)
                  data = compute_dvE(state);
                elseif(mf_display(i)==p_i)
                  data = compute_p_i(state);
                elseif(mf_display(i)==p)
                  data = compute_p(state);
                % Hall terms
                elseif(mf_display(i)==pek)
                  data = compute_pek(state);
                elseif(mf_display(i)==Hal)
                  data = compute_Hal(state);
                elseif(mf_display(i)==Bxu)
                  data = compute_Bxu(state);
                elseif(mf_display(i)==dtJ)
                  data = compute_dtJ(state);
                elseif(mf_display(i)==Eck)
                  % compute the residual of Ohm's law.
                  % There is a lot of redundant computation here
                  data = compute_Bxu(state) ...
                       + compute_Hal(state) ...
                       + compute_pek(state) ...
                       + compute_dtJ(state) ...
                       - compute_E(state,E);
                  data(:,:,3) = sqrt(data(:,:,1).^2 + data(:,:,2).^2 + data(:,:,3).^2);
                elseif(mf_display(i)==Mi_t)
                  data = compute_Mi_t(state);
                  data(:,:,3) = sqrt(data(:,:,1).^2+data(:,:,2).^2+data(:,:,3).^2);
                else
                   error(['support for ' all_names(i) 'not yet implemented']);
                end
              end
              toc;
              %
              % display the data
              disp(['displaying data']);
              tic;
              if(use_subplots)
                subplot(size(mf_display,1),size(mf_display,2),i);
              else
                fig=mf_display(i);
                if(ishandle(fig))
                  set(0,'CurrentFigure',fig);
                else
                  figure(fig);
                end
                clf;
              end
              %
              hold on;
              %
              var_name = char(names(i));
              if(num_components(i)==3)
               % display the third component with color
                plot_scalar(data(:,:,3),xl,yl,flips(i),flips(i),color_ranges(:,i));
                %plot_3d_vector(xl,yl,data(:,:,3),xc,yc,...
                %   data(:,:,1),data(:,:,2),time,var_name,...
                %   flips(i),color_ranges(:,i),vector_scales(i),noDiv(i));
              end
              if(num_components(i)==2 || num_components(i)==3)
                 if(noDiv(i)==1 || noDiv(i)==3)
                   plot_fieldlines(xc,yc,data(:,:,1),data(:,:,2),...
                      params.dx,params.dy);
                 elseif(noDiv(i)==2 || noDiv(i)==3)
                   plot_vector(xc,yc,data(:,:,1),data(:,:,2),flips(i),0);
                 else
                 %   % show nothing for the vector component
                 end
              elseif(num_components(i)==1)
              % display with color
                 plot_scalar(data(:,:,1),xl,yl,... %time,var_name,...
                    flips(i),flips(i),color_ranges(:,i));
              else
                 error([' invalid num_components: ', num2str(num_components(i))]);
              end
              set_axes(domain_scaling);
              set_title(var_name,time);
              hold off;
              toc;
           end
        end
      end
    end
  end
  disp(' ')
end

function params = get_parameters(outputdir)
  parameters_ini = [outputdir '/parameters.ini'];
  params = read_params_from_ini(parameters_ini, ...
    ... % dogpack parameters
    'nout', 'tfinal', 'dtv', 'cflv', 'nv', 'space_order', ...
    'meqn', 'mx', 'my', 'mbc', 'xlow', 'xhigh', 'ylow', 'yhigh', ...
    'mrestart', 'nstart', 'datafmt', ...
    ... % plasma parameters
    'gamma','B_0','BCs','domain_scaling','enforced_symmetry', ...
    'mass_ratio', 'spc_mass_mode', 'ion_mass', 'total_mass', ...
    'temp_ratio','cc','cs_light',...
    'ion_iso_period','elc_iso_period');

  % check parameters
  assert(params.temp_ratio==1);
  assert(params.mass_ratio==1);

  % reset dogpack parameters
  params.xhigh=4*pi*params.domain_scaling;
  params.xlow=-params.xhigh;
  params.yhigh=2*pi*params.domain_scaling;
  params.ylow=-params.yhigh;
  if(bitand(params.enforced_symmetry,1)==1)
    params.xlow=0;
  end
  if(bitand(params.enforced_symmetry,2)==2)
    params.ylow=0;
  end
  % bit 4 means mx and my describe actual computational mesh
  if(bitand(params.enforced_symmetry,5)==1)
    assert(bitand(params.mx,1)==0); % mx is even
    params.mx = params.mx/2;
  end
  if(bitand(params.enforced_symmetry,6)==2)
    assert(bitand(params.my,1)==0); % my is even
    params.my = params.my/2;
  end

  % define derived  application parameters
  %
  params.one_over_epsilon = params.cs_light*params.cs_light;
  params.dx = (params.xhigh - params.xlow)/params.mx;
  params.dy = (params.yhigh - params.ylow)/params.my;
  params.plot_dx = params.dx/params.space_order;
  params.plot_dy = params.dy/params.space_order;
  [ion_mass, elc_mass] = get_species_masses(params);
  params.ion_mass = ion_mass;
  params.elc_mass = elc_mass;
end

function [ion_mass, elc_mass] = get_species_masses(params)
  mass_ratio = params.mass_ratio;
  if(isfield(params,'spc_mass_mode') & params.spc_mass_mode=='ion_mass')
    ion_mass = 1.;
    if(isfield(params,'ion_mass')) ion_mass = params.ion_mass; end;
    elc_mass = ion_mass/mass_ratio;
  else
    total_mass = 1.;
    if(isfield(params,'total_mass')) total_mass = params.total_mass; end;
    elc_mass = total_mass/(mass_ratio+1.);
    ion_mass = mass_ratio*elc_mass;
  end
end

function [data,time]=read_from_file(datafmt,basefilename,mx,my,method1,num_components)
  [data,time]=read_output(datafmt, basefilename, ...
    mx, my, method1, num_components, 1:num_components,1);
end

function [xl,yl,xc,yc]=get_grid2(params)
  [xl,yl]=get_left_grid2(params);
  % centers of cells
  xc=xl(1:end-1,1:end-1)+(0.5*params.plot_dx);
  yc=yl(1:end-1,1:end-1)+(0.5*params.plot_dy);
end

function [xl,yl]=get_left_grid2(params)
  % low/left edges of cells
  xl=params.xlow+(0:params.mx*params.space_order)*(params.plot_dx);
  yl=params.ylow+(0:params.my*params.space_order)*(params.plot_dy);
  % should the user do this?
  [xl,yl]=ndgrid(xl,yl);
end

