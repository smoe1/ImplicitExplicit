function plotout(outputdir_in)

  setarg_outputdir;
  
  all_names = { ...
     ... % state variables
     'rho_i', 'Mi', 'Ni', 'rho_e', 'Me', 'Ne', ...
     'B', 'E', 'Psi', 'phi', ...
     'bogus', ...
      ... % quantities to be constructed
     'sgm', 'dvB', 'dvE', 'u', 'J', 'Ji', 'Je', ...
     'Pi', 'Pe', 'P', 'p_i', 'p_e', 'p', ...
     ... % Ohm's law terms
     'Bxu', 'Hal', 'pek', 'J_t', 'DFJ', 'dtJ', 'Eck', ...
     ... % Ohm's law--related terms
     'Mi_t', 'Me_t', ...
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
  all_num_components(index.Ni      ) = 6;
  all_num_components(index.rho_e   ) = 1;
  all_num_components(index.Me      ) = 3;
  all_num_components(index.Ne      ) = 6;
  all_num_components(index.B       ) = 3;
  all_num_components(index.E       ) = 3;
  all_num_components(index.Psi     ) = 1;
  all_num_components(index.phi     ) = 1;
  %
  all_num_components(index.Bxu) = 3;
  all_num_components(index.Hal) = 3;
  all_num_components(index.pek) = 3; % change to 2?
  all_num_components(index.J_t) = 3;
  all_num_components(index.DFJ) = 3;
  all_num_components(index.dtJ) = 3;
  all_num_components(index.Eck) = 3;
  all_num_components(index.sgm) = 1;
  all_num_components(index.dvB) = 1;
  all_num_components(index.dvE) = 1;
  all_num_components(index.u  ) = 3;
  all_num_components(index.J  ) = 3;
  all_num_components(index.Ji ) = 3;
  all_num_components(index.Je ) = 3;
  all_num_components(index.Pi ) = 6;
  all_num_components(index.Pe ) = 6;
  all_num_components(index.P  ) = 6;
  all_num_components(index.p_i) = 1;
  all_num_components(index.p_e) = 1;
  all_num_components(index.p  ) = 1;
  all_num_components(index.Mi_t)= 3;
  all_num_components(index.Me_t)= 3;
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
  all_color_ranges(:,index.rho_e   ) = [0,0];
  all_color_ranges(:,index.Me      ) = [-0,0];
  all_color_ranges(:,index.Ne      ) = [-0,1];
  all_color_ranges(:,index.B       ) = [-.2,.2];
  all_color_ranges(:,index.E       ) = [-.2,.2];
  all_color_ranges(:,index.Psi     ) = [-.05,.05];
  all_color_ranges(:,index.phi     ) = [-0,0];
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
  all_color_ranges(:,index.Je      ) = [-0,0];
  all_color_ranges(:,index.p       ) = [ 0,1];
  all_color_ranges(:,index.P       ) = [ 0,1];

  all_vector_scales=zeros(1,length(all_names));
  all_vector_scales(:,index.rho_i   ) = 0;
  all_vector_scales(:,index.Mi      ) = all_color_ranges(2,index.Mi);
  all_vector_scales(:,index.Ni      ) = 0;
  all_vector_scales(:,index.rho_e   ) = 0;
  all_vector_scales(:,index.Me      ) = all_color_ranges(2,index.Me);
  all_vector_scales(:,index.Ne      ) = 0;
  all_vector_scales(:,index.B       ) = 1;
  all_vector_scales(:,index.E       ) = all_color_ranges(2,index.E);
  all_vector_scales(:,index.Psi     ) = 0;
  all_vector_scales(:,index.phi     ) = 0;
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
                elseif(mf_display(i)==Je)
                  data = compute_Je(state);
                elseif(mf_display(i)==dvB)
                  data = compute_dvB(state);
                elseif(mf_display(i)==dvE)
                  data = compute_dvE(state);
                elseif(mf_display(i)==Pi)
                  data = compute_Pi(state);
                elseif(mf_display(i)==Pe)
                  data = compute_Pe(state);
                elseif(mf_display(i)==P)
                  data = compute_P(state);
                elseif(mf_display(i)==p_i)
                  data = compute_p_i(state);
                elseif(mf_display(i)==p_e)
                  data = compute_p_e(state);
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
                elseif(mf_display(i)==Me_t)
                  data = compute_Me_t(state);
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
              elseif(num_components(i)==6)
                 plot_eigs(data,xl,yl,flips(i),color_ranges(:,i));
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

  check_parameters(params);

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

function check_parameters(params)
  % check parameters
  %
  % should not be testing for strict inequality for floating point
  if(~(test_equal(params.temp_ratio^2, params.mass_ratio)))
    warndlg(['temperature ratio ' num2str(params.temp_ratio)...
      ' should normally be the square root of the mass ratio '...
      num2str(params.mass_ratio)]);
  end
  if(params.ion_iso_period<=0) assert(params.elc_iso_period<=0);
  elseif(~test_equal(params.ion_iso_period, ...
                 params.elc_iso_period*sqrt(params.mass_ratio)))
    warndlg(['the ratio of the ion_iso_period ' num2str(params.ion_iso_period)...
      ' to the elc_iso_period should normally be' ...
      ' the square root of the mass ratio '...
      num2str(params.mass_ratio)]);
  end
end

% test if two floats are equal
%
function bool = test_equal(num1, num2)
  bool=false;
  eps=1e-13;
  mag=abs(num1)+abs(num2);
  diff=abs(num1-num2);
  if(mag<1)
    if(diff<eps) bool=true; end;
  else
    if(diff/mag < eps) bool = true; end
  end
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

% conversion works by:
%
% state_lc = project_onto_legendre_basis(state_qv,space_order);
% state_qv = sample_state2(state_lc, space_order, 2);
% state_rv = sample_state2(state_lc, space_order, 1);
%
% data structure
%
%
%function qv_to_lc
%end
%
%function lc_to_qv
%end
%
%function to_rv
%end

function data_qv = convert_lc_to_qv(data_lc,space_order)
  data_qv = sample_state2(data_lc, space_order, 2);
end

function data_rv = convert_qv_to_rv(data_qv, space_order)
  data_lc = project_onto_legendre_basis(data_qv, space_order);
  data_rv = sample_state2(data_lc, space_order, 1);
end

function data = compute_sgm(state)
  global space_order state_components rho_i rho_e ion_mass elc_mass;
  sgm_state = state(:,:,:,state_components{rho_i})/ion_mass ...
            - state(:,:,:,state_components{rho_e})/elc_mass;
  data = sample_state2(sgm_state,space_order);
end

function data = compute_u(state)
  global space_order state_components rho_i rho_e Mi Me;
  M_state = state(:,:,:,state_components{Mi}) ...
          + state(:,:,:,state_components{Me});
  rho_state = state(:,:,:,state_components{rho_i}) ...
            + state(:,:,:,state_components{rho_e});
  data = compute_u_val(M_state,rho_state,space_order);
end

function data = compute_J(state)
  global space_order state_components Mi Me
  global ion_mass elc_mass;
  J_state = state(:,:,:,state_components{Mi})/ion_mass ...
          - state(:,:,:,state_components{Me})/elc_mass;
  data = sample_state2(J_state,space_order);
end

function data = compute_Ji(state)
  global space_order state_components Mi
  global ion_mass
  Ji_state = state(:,:,:,state_components{Mi})/ion_mass;
  data = sample_state2(Ji_state,space_order);
end

function data = compute_Je(state)
  global space_order state_components Me
  global elc_mass;
  Je_state = - state(:,:,:,state_components{Me})/elc_mass;
  data = sample_state2(Je_state,space_order);
end

function data = compute_dvB(state)
  global space_order state_components B;
  global dx dy;

  % compute which components have flip symmetry
  %
  global BCs enforced_symmetry;

  % === horizontal symmetries ===
  % mirror boundary conditions (if there is enforced symmetry in the x-axis)
  %  flip_indices: _B2   , _B3

  % === vertical symmetries ===
  %
  % types of boundary conditions
  %
  CONDUCTING_WALL=1;
  PERIODIC=2;
  PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN=3;
  OPEN_BOUNDARY=4;
  %
  % // bottom mirror boundary for enforced symmetry,
  % // top boundary also for periodic and enforced symmetry
  % setBottomTopSymmetryBoundaries
  %     flip_indices: _B1   , _B3
  % // CONDUCTING_WALL
  % setBottomTopConductingWallBoundaries
  %     flip_indices: _B2
  % // PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
  % setBottomTopPeriodicBoundaries
  %     flip_indices: _B2   , _B3
  %

  % magnetic field symmetries
  %
  % x-direction (yes if B1 flipped in horizontal symmetries)
  flipx1 = [1, 1];
  % y-direction (yes if B2 flipped in vertical symmetries)
  if(BCs==PERIODIC)
    flipy2 = [1, 1];
  else
    flipy2 = [-1,-1];
  end
  if(bitand(enforced_symmetry,2)==2)
    flipy2(1) = 1;
  end

  B_state = state(:,:,:,state_components{B});
  div_B_state = div_stateVector(B_state,space_order,dx,dy,flipx1,flipy2);
  data = sample_state2(div_B_state,space_order);
end

function data = compute_dvE(state)
  global space_order state_components E;
  global dx dy;

  % compute which components have flip symmetry
  %
  global BCs enforced_symmetry;

  % === horizontal symmetries ===
  % mirror boundary conditions (if there is enforced symmetry in the x-axis)
  %  flip_indices: _E1

  % === vertical symmetries ===
  %
  % types of boundary conditions
  %
  CONDUCTING_WALL=1;
  PERIODIC=2;
  PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN=3;
  OPEN_BOUNDARY=4;
  %
  % // bottom mirror boundary for enforced symmetry,
  % // top boundary also for periodic and enforced symmetry
  % setBottomTopSymmetryBoundaries
  %     flip_indices: _E2
  % // CONDUCTING_WALL
  % setBottomTopConductingWallBoundaries
  %     flip_indices: _E1   , _E3
  % // PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
  % setBottomTopPeriodicBoundaries
  %     flip_indices: _E2   , _E3
  %

  % electric field symmetries
  %
  % x-direction (yes if E1 flipped in horizontal symmetries)
  flipx1 = [-1, -1];
  % y-direction (yes if E2 flipped in vertical symmetries)
  if(BCs==CONDUCTING_WALL)
    flipy2 = [1, 1];
  else
    flipy2 = [-1,-1];
  end
  if(bitand(enforced_symmetry,2)==2)
    flipy2(1) = -1;
  end

  E_state = state(:,:,:,state_components{E});
  div_E_state = div_stateVector(E_state,space_order,dx,dy,flipx1,flipy2);
  data = sample_state2(div_E_state,space_order);
end

function data = compute_Pi(state)
  global space_order state_components rho_i Mi Ni;
  Mi_state=state(:,:,:,state_components{Mi});
  rho_i_state = state(:,:,:,state_components{rho_i});
  KEi_state = compute_KE_state(Mi_state,rho_i_state,space_order);
  Ni_state=state(:,:,:,state_components{Ni});
  Pi_state=Ni_state-KEi_state;
  data = sample_state2(Pi_state,space_order);
end

function data = compute_Pe(state)
  global space_order state_components rho_e Me Ne;
  Me_state=state(:,:,:,state_components{Me});
  rho_e_state = state(:,:,:,state_components{rho_e});
  KEe_state = compute_KE_state(Me_state,rho_e_state,space_order);
  Ne_state=state(:,:,:,state_components{Ne});
  Pe_state=Ne_state-KEe_state;
  data = sample_state2(Pe_state,space_order);
end

function data = compute_P(state)
  %
  % ions
  %
  global space_order state_components rho_i Mi Ni;
  Mi_state=state(:,:,:,state_components{Mi});
  rho_i_state = state(:,:,:,state_components{rho_i});
  KEi_state = compute_KE_state(Mi_state,rho_i_state,space_order);
  Ni_state=state(:,:,:,state_components{Ni});
  Pi_state=Ni_state-KEi_state;
  %
  % elcs
  %
  global space_order state_components rho_e Me Ne;
  Me_state=state(:,:,:,state_components{Me});
  rho_e_state = state(:,:,:,state_components{rho_e});
  KEe_state = compute_KE_state(Me_state,rho_e_state,space_order);
  Ne_state=state(:,:,:,state_components{Ne});
  Pe_state=Ne_state-KEe_state;
  %
  P_state = Pi_state + Pe_state;
  data = sample_state2(P_state,space_order);
end

% P11 = S1
% P12 = S2
% P13 = S3
% P22 = S4
% P23 = S5
% P33 = S6
function data = compute_E(state,E)
  global space_order state_components;
  components=state_components{E};
  data=sample_state2(state(:,:,:,components),space_order);
end

function data = compute_p_i(state)
  data = compute_Pi(state);
  data = (data(:,:,1) + data(:,:,4) + data(:,:,6))/3.;
end

function data = compute_p_e(state)
  data = compute_Pe(state);
  data = (data(:,:,1) + data(:,:,4) + data(:,:,6))/3.;
end

function data = compute_p(state)
  data = compute_P(state);
  data = (data(:,:,1) + data(:,:,4) + data(:,:,6))/3.;
end

function data = compute_pek(state)
  global space_order state_components rho_i rho_e Mi Me Ni Ne;
  global ion_mass elc_mass;
  global dx dy;
  % ions
  %
  Mi_state=state(:,:,:,state_components{Mi});
  rho_i_state = state(:,:,:,state_components{rho_i});
  KEi_state = compute_KE_state(Mi_state,rho_i_state,space_order);
  Ni_state=state(:,:,:,state_components{Ni});
  Pi_state=Ni_state-KEi_state;
  %
  % elcs
  %
  Me_state=state(:,:,:,state_components{Me});
  rho_e_state = state(:,:,:,state_components{rho_e});
  KEe_state = compute_KE_state(Me_state,rho_e_state,space_order);
  Ne_state=state(:,:,:,state_components{Ne});
  Pe_state=Ne_state-KEe_state;
  
  % electrokinetic pressure
  %Pek_state=Pi_state/ion_mass-Pe_state/elc_mass;
  Pek_state=Pi_state*elc_mass-Pe_state*ion_mass;
  % take the divergence of the electrokinetic pressure
  divPek_state=div_stateTensor(Pek_state,space_order,dx,dy);
  divPek_qv = sample_state2(divPek_state,space_order,2);
  rho_state = state(:,:,:,state_components{rho_i}) ...
            + state(:,:,:,state_components{rho_e});
  rho_qv = sample_state2(rho_state,space_order,2);
  data_qv = tupleOverScalar(divPek_qv,rho_qv);
  data_lc = project_onto_legendre_basis(data_qv,space_order);
  data = sample_state2(data_lc, space_order, 1);
end

function data = compute_Hal(state)
  global space_order state_components rho_i rho_e Mi Me B;
  global ion_mass elc_mass;

  J_state = state(:,:,:,state_components{Mi})/ion_mass ...
          - state(:,:,:,state_components{Me})/elc_mass;
  J_val = sample_state2(J_state,space_order);
  B_state = state(:,:,:,state_components{B});
  B_val = sample_state2(B_state,space_order);
  JxB_vals = crossProduct(J_val,B_val);
  rho_state = state(:,:,:,state_components{rho_i}) ...
            + state(:,:,:,state_components{rho_e});
  rho_val = sample_state2(rho_state,space_order);
  data = tupleTimesScalar(JxB_vals,(ion_mass-elc_mass)./rho_val);
end

function data = compute_Bxu(state)
  global space_order state_components rho_i rho_e Mi Me B;
  M_state = state(:,:,:,state_components{Mi}) ...
          + state(:,:,:,state_components{Me});
  rho_state = state(:,:,:,state_components{rho_i}) ...
            + state(:,:,:,state_components{rho_e});
  uval = compute_u_val(M_state,rho_state,space_order);
  B_state = state(:,:,:,state_components{B});
  B_val = sample_state2(B_state,space_order);
  data = crossProduct(B_val,uval);
end

function data = compute_dtJ(state)
  global space_order state_components;
  global rho_i rho_e Mi Me Ni Ne B E;
  global ion_mass elc_mass;
  global dx dy;
  % Suffixes indicate one of three ways of representing states:
  %
  % qv: quadrature values.
  % lc: legendre coefficients
  % rg: regular grid samples
  %
  % high-order numerical differentation needs the lc representation:
  %   div_state_lc = div_stateTensor(state_lc,space_order),
  % arithmetic operation needs the qv or rg representation, and
  % for output we use the rg representation.
  %
  % conversions supported are qv <-> lc -> rg, i.e.,:
  %   qv -> lc
  %   lc -> qv
  %   lc -> rg
  %
  % Essentially we need to stick with the lc/qv representation
  % until we are done with differentiation, at which point we
  % can move to rg representation.  So if we want to calculate
  % the *curl* of the Ohm's law terms, we will have to delay
  % moving to the rg representation.
  %
  % conversion works by:
  %
  % state_lc = project_onto_legendre_basis(state_qv,space_order);
  % state_qv = sample_state2(state_lc, space_order, 2);
  % state_rv = sample_state2(state_lc, space_order, 1);

  % the multiplier
  %
  rho_i_lc = state(:,:,:,state_components{rho_i});
  rho_e_lc = state(:,:,:,state_components{rho_e});
  rho_i_qv = sample_state2(rho_i_lc,space_order,2);
  rho_e_qv = sample_state2(rho_e_lc,space_order,2);
  rho_qv = rho_i_qv + rho_e_qv;
  mass_prod = ion_mass*elc_mass;
  mass_prod_ovr_rho = mass_prod./rho_qv;

  % for the source term
  %
  E_lc = state(:,:,:,state_components{E});
  E_qv = sample_state2(E_lc,space_order,2);
  B_lc = state(:,:,:,state_components{B});
  B_qv = sample_state2(B_lc,space_order,2);

  % ions
  %
  % flux term
  %
  Ni_state=state(:,:,:,state_components{Ni});
  Mi_flux_lc = div_stateTensor(Ni_state,space_order,dx,dy);
  Mi_flux_qv = sample_state2(Mi_flux_lc, space_order, 2);
  %
  % source term
  %
  sgm_i_qv = rho_i_qv/ion_mass;
  E_tms_sgm_i_qv = tupleTimesScalar(E_qv,sgm_i_qv);
  Mi_lc = state(:,:,:,state_components{Mi});
  Ji_qv = sample_state2(Mi_lc/ion_mass, space_order, 2);
  J_i_cross_B_qv = crossProduct(Ji_qv,B_qv);
  source_i_qv = E_tms_sgm_i_qv + J_i_cross_B_qv;
  % total forcing
  Mi_t_qv = source_i_qv - Mi_flux_qv;

  % elcs
  %
  % flux term
  %
  Ne_state=state(:,:,:,state_components{Ne});
  Me_flux_lc = div_stateTensor(Ne_state,space_order,dx,dy);
  Me_flux_qv = sample_state2(Me_flux_lc, space_order, 2);
  %
  % source term
  %
  sgm_e_qv = -rho_e_qv/elc_mass;
  E_tms_sgm_e_qv = tupleTimesScalar(E_qv,sgm_e_qv);
  Me_lc = state(:,:,:,state_components{Me});
  Je_qv = sample_state2(-Me_lc/elc_mass, space_order, 2);
  J_e_cross_B_qv = crossProduct(Je_qv,B_qv);
  source_e_qv = E_tms_sgm_e_qv + J_e_cross_B_qv;
  % total forcing
  Me_t_qv = source_e_qv - Me_flux_qv;

  J_t_qv = Mi_t_qv/ion_mass - Me_t_qv/elc_mass;

  % DivJflux_lc
  %
  M_state = state(:,:,:,state_components{Mi}) ...
          + state(:,:,:,state_components{Me});
  rho_state = state(:,:,:,state_components{rho_i}) ...
            + state(:,:,:,state_components{rho_e});
  u_qv = compute_u_val(M_state,rho_state,space_order,2);
  J_state = state(:,:,:,state_components{Mi})/ion_mass ...
          - state(:,:,:,state_components{Me})/elc_mass;
  J_qv = sample_state2(J_state,space_order,2);
  uJ_plus_Ju_qv = twiceSymmetricTensorProduct(u_qv, J_qv);
  JJ_qv = selfSymmetricTensorProduct(J_qv);
  JJ_term_qv = tupleTimesScalar(JJ_qv,(ion_mass-elc_mass)./rho_qv);
  Jflux_qv = uJ_plus_Ju_qv - JJ_term_qv;
  Jflux_lc = project_onto_legendre_basis(Jflux_qv,space_order);
  DivJflux_lc = div_stateTensor(Jflux_lc,space_order,dx,dy);
  DivJflux_qv = sample_state2(DivJflux_lc, space_order, 2);

  J_t_term_qv = tupleTimesScalar(J_t_qv,mass_prod_ovr_rho);
  % why do I fail to get agreement when I include DivJflux_qv?
  intertial_term_qv = J_t_term_qv; % + DivJflux_qv;

  % this step would be unnecessary if we had worked with
  % rg values, but I wanted to retain the possibility of converting
  % the results to lc representation so we could take the curl
  data = convert_qv_to_rv(intertial_term_qv, space_order);
end

function data = compute_Mi_t(state)
  global space_order state_components;
  global rho_i Mi Ni B E;
  global ion_mass;
  global dx dy;

  % the multiplier
  %
  rho_i_lc = state(:,:,:,state_components{rho_i});
  rho_i_qv = sample_state2(rho_i_lc,space_order,2);

  % for the source term
  %
  E_lc = state(:,:,:,state_components{E});
  E_qv = sample_state2(E_lc,space_order,2);
  B_lc = state(:,:,:,state_components{B});
  B_qv = sample_state2(B_lc,space_order,2);

  % flux term
  %
  Ni_state=state(:,:,:,state_components{Ni});
  Mi_flux_lc = div_stateTensor(Ni_state,space_order,dx,dy);
  Mi_flux_qv = sample_state2(Mi_flux_lc, space_order, 2);
  %
  % source term
  %
  sgm_i_qv = rho_i_qv/ion_mass;
  E_tms_sgm_i_qv = tupleTimesScalar(E_qv,sgm_i_qv);
  Mi_lc = state(:,:,:,state_components{Mi});
  Ji_qv = sample_state2(Mi_lc/ion_mass, space_order, 2);
  J_i_cross_B_qv = crossProduct(Ji_qv,B_qv);
  source_i_qv = E_tms_sgm_i_qv + J_i_cross_B_qv;
  % total forcing
  Mi_t_qv = source_i_qv - Mi_flux_qv;
  % this step would be unnecessary if we had worked with rg values
  data = convert_qv_to_rv(Mi_t_qv, space_order);
end

function data = compute_Me_t(state)
  global space_order state_components;
  global rho_e Me Ne B E;
  global elc_mass;
  global dx dy;

  % the multiplier
  %
  rho_e_lc = state(:,:,:,state_components{rho_e});
  rho_e_qv = sample_state2(rho_e_lc,space_order,2);

  % for the source term
  %
  E_lc = state(:,:,:,state_components{E});
  E_qv = sample_state2(E_lc,space_order,2);
  B_lc = state(:,:,:,state_components{B});
  B_qv = sample_state2(B_lc,space_order,2);

  % elcs
  %
  % flux term
  %
  Ne_state=state(:,:,:,state_components{Ne});
  Me_flux_lc = div_stateTensor(Ne_state,space_order,dx,dy);
  Me_flux_qv = sample_state2(Me_flux_lc, space_order, 2);
  %
  % source term
  %
  sgm_e_qv = -rho_e_qv/elc_mass;
  E_tms_sgm_e_qv = tupleTimesScalar(E_qv,sgm_e_qv);
  Me_lc = state(:,:,:,state_components{Me});
  Je_qv = sample_state2(-Me_lc/elc_mass, space_order, 2);
  J_e_cross_B_qv = crossProduct(Je_qv,B_qv);
  source_e_qv = E_tms_sgm_e_qv + J_e_cross_B_qv;
  % total forcing
  Me_t_qv = source_e_qv - Me_flux_qv;
  % this step would be unnecessary if we had worked with rg values
  data = convert_qv_to_rv(Me_t_qv, space_order);
end

function product = tupleTimesScalar(tuple,scalar)
  product = zeros(size(tuple));
  for idx=1:size(tuple,3)
    product(:,:,idx) = tuple(:,:,idx).*scalar;
  end
end

function product = tupleOverScalar(tuple,scalar)
  product = zeros(size(tuple));
  for idx=1:size(tuple,3)
    product(:,:,idx) = tuple(:,:,idx)./scalar;
  end
end

function out = crossProduct(A,B)
  assert(all(size(A)==size(B)));
  assert(size(size(A),2)==3);
  assert(size(A,3)==3);
  out = zeros(size(A));
  out(:,:,1) = A(:,:,2).*B(:,:,3) - A(:,:,3).*B(:,:,2);
  out(:,:,2) = A(:,:,3).*B(:,:,1) - A(:,:,1).*B(:,:,3);
  out(:,:,3) = A(:,:,1).*B(:,:,2) - A(:,:,2).*B(:,:,1);
end

function uJ_plus_Ju_qv = twiceSymmetricTensorProduct(u_qv, J_qv)
  assert(size(J_qv,3)==3);
  assert(all(size(u_qv)==size(J_qv)));
  dims=size(J_qv);
  tensorDims=[dims(1:2),6];
  % T11 = S1
  % T12 = S2
  % T13 = S3
  % T22 = S4
  % T23 = S5
  % T33 = S6
  uJ_plus_Ju_qv = zeros(tensorDims);
  uJ_plus_Ju_qv(:,:,1) = 2*u_qv(:,:,1).*J_qv(:,:,1);
  uJ_plus_Ju_qv(:,:,2) = u_qv(:,:,1).*J_qv(:,:,2) + u_qv(:,:,2).*J_qv(:,:,1);
  uJ_plus_Ju_qv(:,:,3) = u_qv(:,:,1).*J_qv(:,:,3) + u_qv(:,:,3).*J_qv(:,:,1);
  uJ_plus_Ju_qv(:,:,4) = 2*u_qv(:,:,2).*J_qv(:,:,2);
  uJ_plus_Ju_qv(:,:,5) = u_qv(:,:,2).*J_qv(:,:,3) + u_qv(:,:,3).*J_qv(:,:,2);
  uJ_plus_Ju_qv(:,:,6) = 2*u_qv(:,:,3).*J_qv(:,:,3);
end

function JJ_qv = selfSymmetricTensorProduct(J_qv)
  assert(size(J_qv,3)==3);
  dims=size(J_qv);
  tensorDims=[dims(1:2),6];
  JJ_qv = zeros(tensorDims);
  JJ_qv(:,:,1) = J_qv(:,:,1).*J_qv(:,:,1);
  JJ_qv(:,:,2) = J_qv(:,:,1).*J_qv(:,:,2);
  JJ_qv(:,:,3) = J_qv(:,:,1).*J_qv(:,:,3);
  JJ_qv(:,:,4) = J_qv(:,:,2).*J_qv(:,:,2);
  JJ_qv(:,:,5) = J_qv(:,:,2).*J_qv(:,:,3);
  JJ_qv(:,:,6) = J_qv(:,:,3).*J_qv(:,:,3);
end

function div_state=div_stateVector(state,space_order,dx,dy,flipx1,flipy2)

  % compute which components have flip symmetry
  %
  global BCs enforced_symmetry;

  % === horizontal symmetries ===
  % mirror boundary conditions (if there is enforced symmetry in the x-axis)
  %  flip_indices: _M1_i , _M1_e , _B2   , _B3   , _E1   , _psi

  % === vertical symmetries ===
  % // bottom mirror boundary for enforced symmetry,
  % // top boundary also for periodic and enforced symmetry
  % setBottomTopSymmetryBoundaries
  %     flip_indices: _M2_i , _M2_e , _B1   , _B3   , _E2
  % // CONDUCTING_WALL
  % setBottomTopConductingWallBoundaries
  %     flip_indices: _M2_i , _M2_e , _B2   , _E1   , _E3
  % // PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
  % setBottomTopPeriodicBoundaries
  %     flip_indices: _M2_i , _M3_i , _M2_e , _M3_e , _B2   , _B3   , _E2   , _E3
  %

  % default: momentum symmetries
  if(nargin<5)
    flipx1 = [1, 1];
    flipy2 = [-1,-1];
  end

  dims=size(state);
  assert(dims(4)==2||dims(4)==3);
  div_state=zeros([dims(1:3),1]);
  div_state(:,:,:,1) = dx_state(state(:,:,:,1),space_order,dx,flipx1) ...
                     + dy_state(state(:,:,:,2),space_order,dy,flipy2);
end

% Take the divergence of a symmetric tensor P = S whose components
% are in the order:
% P11 = S1
% P12 = S2
% P13 = S3
% P22 = S4
% P23 = S5
% P33 = S6
%
% divP_1 = dx_P11 + dy_P21 = dx_S1 + dy_S2
% divP_2 = dx_P12 + dy_P22 = dx_S2 + dy_S4
% divP_3 = dx_P13 + dy_P23 = dx_S3 + dy_S5
%
% (the sixth component is not used)
function div_state=div_stateTensor(state,space_order,dx,dy)
  dims=size(state);
  assert(dims(4)==6);
  div_state=zeros([dims(1:3),3]);

  % compute which components have flip symmetry
  %
  global BCs enforced_symmetry;

  % === horizontal symmetries ===
  %
  % populate_ghost_cells
  %
  % is there enforced symmetry in the x-axis?
  %if(bitand(enforced_symmetry,1)==1)
  %  % mirror boundary conditions
  %  flip_indices: _M1_i , _N12_i, _N13_i,

  % === vertical symmetries ===
  %
  % types of boundary conditions
  %
  CONDUCTING_WALL=1;
  PERIODIC=2;
  PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN=3;
  OPEN_BOUNDARY=4;
  %
  % // bottom mirror boundary for enforced symmetry,
  % // top boundary also for periodic and enforced symmetry
  % setBottomTopSymmetryBoundaries
  %     flip_indices: _M2_i , _N12_i, _N23_i,
  % // CONDUCTING_WALL
  % setBottomTopConductingWallBoundaries
  %     flip_indices: _M2_i , _N12_i, _N23_i,
  % // PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
  % setBottomTopPeriodicBoundaries
  %     flip_indices: _M2_i , _M3_i , _N12_i, _N13_i,
  %
  flipy23 = [-1,-1];
  if(BCs==PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN)
    flip23(2)=1;
    if(bitand(enforced_symmetry,2)==0)
      flip23(1) = 1;
    end
  end

  div_state(:,:,:,1) = dx_state(state(:,:,:,1),space_order,dx,[ 1, 1]) ... % N11
                     + dy_state(state(:,:,:,2),space_order,dy,[-1,-1]);    % N21
  div_state(:,:,:,2) = dx_state(state(:,:,:,2),space_order,dx,[-1,-1]) ... % N12
                     + dy_state(state(:,:,:,4),space_order,dy,[ 1, 1]);    % N22
  div_state(:,:,:,3) = dx_state(state(:,:,:,3),space_order,dx,[-1,-1]) ... % N13
                     + dy_state(state(:,:,:,5),space_order,dy,flipy23);    % N23
end

% handle mirroring and wrap-around without actually padding the array.
function out = dx_state(state, space_order,dx,flip)
  assert(all(flip.*flip==[1 1]));
  assert(space_order >= 1);
  assert(space_order <= 3);
  dims=size(state);
  assert(dims(3)==get_kmax(space_order));
  out = zeros(dims);
  paddims = dims;
  paddims(1) = 1;

  % determine ghost cell contents
  % based on enforced symmetries and the periodic boundary conditions
  %
  global enforced_symmetry;
  leftpad = zeros(paddims);
  rghtpad = zeros(paddims);
  if(bitand(enforced_symmetry,1)==1) % enforced symmetry in X axis?
    reversed_components = [1, -1,1, -1,1,1]; % 1st, 2nd, 3rd order components
    left_sgns = reversed_components*flip(1);
    rght_sgns = reversed_components*flip(2);
    for i=1:get_kmax(space_order)
      leftpad(1,:,i,:) = left_sgns(i)*state(  1,:,i,:);
      rghtpad(1,:,i,:) = rght_sgns(i)*state(end,:,i,:);
    end
  else % periodic
    % for the periodic case simply copy wrapped array elements
    leftpad(1,:,:,:) = state(end,:,:,:);
    rghtpad(1,:,:,:) = state(  1,:,:,:);
  end

  % take centered differences of all components
  out(2:end-1,:,:,:) = (state(3:end,:,:,:)-state(1:end-2,:,:,:))/(2.*dx);
  out(1,:,:,:) = (state(2,:,:,:)-leftpad)/(2.*dx);
  out(end,:,:,:) = (rghtpad - state(end-1,:,:,:))/(2.*dx);
  if(space_order>=3)
    % correct first component for third-order accuracy
    out(2:end-1,:,1,:) = out(2:end-1,:,1,:) ...
      - (state(3:end,:,5,:)-state(1:end-2,:,5,:))*(sqrt(5)/dx);
    out(1,:,1,:) = out(1,:,1,:) ...
      - (state(2,:,5,:)-leftpad(1,:,5,:))*(sqrt(5)/dx);
    out(end,:,1,:) = out(1,:,1,:) ...
      - (rghtpad(1,:,5,:) - state(end-1,:,5,:))*(sqrt(5)/dx);
  end
end

function out = dy_state(state, space_order,dy,flip)
  assert(all(flip.*flip==[1 1]));
  assert(space_order >= 1);
  assert(space_order <= 3);
  dims=size(state);
  assert(dims(3)==get_kmax(space_order));
  out = zeros(dims);
  paddims = dims;
  paddims(2) = 1;
  %
  % determine ghost cell sources and negations
  % based on boundary conditions and enforced symmetries
  global BCs enforced_symmetry;
  %
  % types of boundary conditions
  %
  CONDUCTING_WALL=1;
  PERIODIC=2;
  PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN=3;
  OPEN_BOUNDARY=4;
  %
  % determine ghost cell contents
  % based on enforced symmetries and the periodic boundary conditions
  %
  low_pad = zeros(paddims);
  highpad = zeros(paddims);
  %
  reversed_components = [1, 1,-1, -1,1,1];
  low__sgns = reversed_components*flip(1);
  high_sgns = reversed_components*flip(2);
  %
  if(bitand(enforced_symmetry,2)==2) % enforced symmetry in Y axis?
    % low end obeys mirror symmetry
    disp('enforcing vertical mirror symmetry');
    for i=1:get_kmax(space_order)
      low_pad(:,1,i,:) = low__sgns(i)*state(:,    1,i,:);
    end
    if(BCs==OPEN_BOUNDARY)
      % copy all components from neighbor without negation
      highpad(:,1,:,:) = state(:,end,:,:);
    else
      for i=1:get_kmax(space_order)
        highpad(:,1,i,:) = high_sgns(i)*state(:,end,i,:);
      end
    end
  else
    if(BCs==PERIODIC)
      % for the periodic case simply copy wrapped array elements
      low_pad(:,1,:,:) = state(:,end,:,:);
      highpad(:,1,:,:) = state(:,  1,:,:);
    elseif(BCs==OPEN_BOUNDARY)
      % copy all components from neighbor without negation
      low_pad(:,1,:,:) = state(:,  1,:,:);
      highpad(:,1,:,:) = state(:,end,:,:);
    else
      for i=1:get_kmax(space_order)
        low_pad(:,1,i,:) = low__sgns(i)*state(:,   1,i,:);
        highpad(:,1,i,:) = high_sgns(i)*state(:,end,i,:);
      end
    end
  end

  % take centered differences of all components
  %
  out(:,2:end-1,:,:) = (state(:,3:end,:,:)-state(:,1:end-2,:,:))/(2.*dy);
  out(:,1,:,:) = (state(:,2,:,:)-low_pad)/(2.*dy);
  out(:,end,:,:) = (highpad - state(:,end-1,:,:))/(2.*dy);
  if(space_order>=3)
    % correct first component for third-order accuracy
    out(:,2:end-1,1,:) = out(:,2:end-1,1,:) ...
      - (state(:,3:end,6,:)-state(:,1:end-2,6,:))*(sqrt(5)/dy);
    out(:,1,1,:) = out(:,1,1,:) ...
      - (state(:,2,6,:)-low_pad(:,1,6,:))*(sqrt(5)/dy);
    out(:,end,1,:) = out(:,end,1,:) ...
      - (highpad(:,1,6,:) - state(:,end-1,6,:))*(sqrt(5)/dy);
  end
end

% Calculate the kinetic energy tensor KE = S whose components
% are in the order:
% KE11 = S1
% KE12 = S2
% KE13 = S3
% KE22 = S4
% KE23 = S5
% KE33 = S6
%
function KE_state = compute_KE_state(M_state,rho_state,space_order)
  % (1) sample state variables at Gaussian quadrature points
  M_quad_vals=sample_state2(M_state,space_order,2);
  rho_quad_vals=sample_state2(rho_state,space_order,2);
  % (2) calculate values at points
  dims=size(rho_quad_vals);
  M2_quad_vals = zeros([dims(1:2), 6]);
  M2_quad_vals(:,:,1) = M_quad_vals(:,:,1).*M_quad_vals(:,:,1);
  M2_quad_vals(:,:,2) = M_quad_vals(:,:,1).*M_quad_vals(:,:,2);
  M2_quad_vals(:,:,3) = M_quad_vals(:,:,1).*M_quad_vals(:,:,3);
  M2_quad_vals(:,:,4) = M_quad_vals(:,:,2).*M_quad_vals(:,:,2);
  M2_quad_vals(:,:,5) = M_quad_vals(:,:,2).*M_quad_vals(:,:,3);
  M2_quad_vals(:,:,6) = M_quad_vals(:,:,3).*M_quad_vals(:,:,3);
  KE_quad_vals = zeros(size(M2_quad_vals));
  for idx=1:6
    KE_quad_vals(:,:,idx) = M2_quad_vals(:,:,idx)./rho_quad_vals;
  end
  % (3) project onto the Legendre basis
  % using Gaussian quadrature for the integration.
  KE_state = project_onto_legendre_basis(KE_quad_vals,space_order);
end

function state = project_onto_legendre_basis(quad_vals,space_order)
  wght_2d = get_wght_2d(space_order);
  mx=size(quad_vals,1)/space_order;
  my=size(quad_vals,2)/space_order;
  kmax=get_kmax(space_order);
  meqns=size(quad_vals,3);
  phi = sample_basis_functions2(space_order, 2);
  for k_idx=1:kmax
    accumulator = zeros(mx,my,meqns);
    for n1=1:space_order
    for n2=1:space_order
      % I guess that the legendre functions are orthonormal.
      % so to get the coefficients of the legendre function expansion
      % we just integrate the state values multiplied by the legendre 
      % function values.  (We integrate by sampling at Gaussian
      % quadrature points and taking a linear combination with
      % Gaussian weights.)
      accumulator = accumulator + wght_2d(n1,n2)*phi(n1,n2,k_idx) ...
        *quad_vals(n1:space_order:end,n2:space_order:end,:);
    end
    end
    state(:,:,k_idx,:) = reshape(accumulator,[mx,my,1,meqns])/4.0;
  end
end

function wght_2d = get_wght_2d(space_order)
  wght_1d = get_wght_1d(space_order);
  wght_2d = transpose(wght_1d) * wght_1d;
end

function wght_1d = get_wght_1d(space_order)
    if(space_order==1)
      wght_1d = 2.0;
    elseif(space_order==2)
      wght_1d = [1., 1.];
    elseif(space_order==3)
      wght_1d = [5., 8., 5.]/9.0;
    elseif(space_order==4)
      w1 = (18.0 - sqrt(3)*sqrt(10))/36.0;
      w2 = (18.0 + sqrt(3)*sqrt(10))/36.0;
      wght_1d = [w1, w2, w2, w1]; 
    elseif(space_order==5)
      w1 = (322.0 - 13.0*sqrt(7)*sqrt(10))/900.0;
      w2 = (322.0 + 13.0*sqrt(7)*sqrt(10))/900.0;
      w3 = 128.0/225.0;
      wght_1d = [w1, w2, w3, w2, w1]; 
    else
      error(['invalid space_order: ' num2str(space_order)]);
    end
end

function u_val = compute_u_val(M_state,rho_state,space_order,point_type);
  if(nargin<4)
    point_type=1;
  end
  M_val = sample_state2(M_state,space_order,point_type);
  rho_val = sample_state2(rho_state,space_order,point_type);
  u_val = zeros(size(M_val));
  for j=1:3 
    u_val(:,:,j) = M_val(:,:,j)./rho_val;
  end
end

function [data,time]=read_from_file(datafmt,basefilename,mx,my,space_order,num_components)

  [data,time]=read_output(datafmt, basefilename, ...
    mx, my, space_order, num_components, 1:num_components,1);
end

% isotropization results in black and white
function plot_eigs(data,xl,yl,flips,color_range)
    % assume data represents components of symmetric tensor
    % and compute eigenvalues
    %
    [eig1, eig2, eig3] = get_eigs(data);
    all_eigs = zeros(size(eig1,2),size(eig1,1),3);
    all_eigs(:,:,1) = eig1';
    all_eigs(:,:,2) = eig2';
    all_eigs(:,:,3) = eig3';
    % need to rescale.
    % ignoring color_ranges here:
    if(color_range(2)-color_range(1)==0)
      maxval = max(max(max(all_eigs)));
      minval = 0;
    else
      minval=color_range(1);
      maxval=color_range(2);
    end
    % minval becomes 0, maxval becomes 1
    rgb_values = (all_eigs-minval)/(maxval-minval);
    % cap values at 1 and floor them at 0
    rgb_values = max(rgb_values,0);
    rgb_values = min(rgb_values,1);
    xl1d = xl(:,1);
    yl1d = yl(1,:);
    colormap('gray');
    colorbar;
    caxis([minval, maxval]);
    image( xl1d, yl1d,rgb_values);
    image(-xl1d, yl1d,rgb_values);
    image( xl1d,-yl1d,rgb_values);
    image(-xl1d,-yl1d,rgb_values);
end

function [eig1, eig2, eig3] = get_eigs(data)
    P11 = data(:,:,1);
    P12 = data(:,:,2);
    P13 = data(:,:,3);
    P22 = data(:,:,4);
    P23 = data(:,:,5);
    P33 = data(:,:,6);
    %
    % (This uses the formula for the roots of a cubic)
    a1=-P11-P22-P33;
    a2=P11.*P22+P11.*P33+P22.*P33-P12.^2-P13.^2-P23.^2;
    a3=P11.*P23.^2+P22.*P13.^2+P33.*P12.^2 ...
      -2*P12.*P13.*P23-P11.*P22.*P33;
    Q=(3*a2-a1.^2)/9;
    R=(9*a1.*a2-27*a3-2*a1.^3)/54;
    discriminant=Q.^3+R.^2;
    % discriminant should be negative
    %find(discriminant>0)
    root_disc = sqrt(discriminant); % should be imaginary
    S=(R+root_disc).^(1/3);
    % could just let T be the complex conjugate of S.
    T=(R-root_disc).^(1/3);
    % should be real
    S_plus_T=real(S+T);
    imag_S_minus_T = imag(S-T);
    eig3=S_plus_T-a1/3;
    partA = -S_plus_T/2-a1/3;
    partB = (sqrt(3)/2)*imag_S_minus_T;
    eig2= partA+partB;
    eig1= partA-partB;
end
