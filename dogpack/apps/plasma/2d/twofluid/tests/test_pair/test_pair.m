% script I used to verify that symmetric pair code
% exactly agrees with non-pair code.
%
% First execute the following commands:
%
%   DOGPACK=~/dogpack
%   TWOFLUID=$(DOGPACK)/apps/plasma/2d/twofluid
%   TEST_PAIR=$(TWOFLUID)/tests/test_pair
%
%   # overwrite parameters.ini files and makefiles
%   # with agreeing symmetric pair plasma versions
%
%   cp $(TEST_PAIR)/p05_parameters.ini $(TWOFLUID)/p05/parameters.ini
%   cp $(TEST_PAIR)/g05_parameters.ini $(TWOFLUID)/g05/parameters.ini
%   cp $(TEST_PAIR)/p10_parameters.ini $(TWOFLUID)/p10/parameters.ini
%   cp $(TEST_PAIR)/g10_parameters.ini $(TWOFLUID)/g10/parameters.ini
%
%   cp $(TEST_PAIR)/p05_Makefile $(TWOFLUID)/p05/Makefile
%   cp $(TEST_PAIR)/g05_Makefile $(TWOFLUID)/g05/Makefile
%   cp $(TEST_PAIR)/p10_Makefile $(TWOFLUID)/p10/Makefile
%   cp $(TEST_PAIR)/g10_Makefile $(TWOFLUID)/g10/Makefile
%
%   # recompile
%
%   cd $(TWOFLUID)/p05; make cleanallo
%   cd $(TWOFLUID)/g05; make cleanallo
%   cd $(TWOFLUID)/p10; make cleanallo
%   cd $(TWOFLUID)/g10; make cleanallo
%
%   cd $(TWOFLUID)/p05; make -j6
%   cd $(TWOFLUID)/g05; make -j6
%   cd $(TWOFLUID)/p10; make -j6
%   cd $(TWOFLUID)/g10; make -j6
%
%   cd $(TWOFLUID)/p05; ./dog.exe
%   cd $(TWOFLUID)/g05; ./dog.exe
%   cd $(TWOFLUID)/p10; ./dog.exe
%   cd $(TWOFLUID)/g10; ./dog.exe
%
function compare()
  vars = { ...
    'Psi'  , ...
    'E3'   , ...
    'B2'   , ...
    'B1'   , ...
    ... %'nrg_i'  , ...
    'Mi3' , ...
    'Mi2' , ...
    'Mi1' , ...
    'rho_i'};
  models={'p10','g10'};
  for i=1:numel(vars)
    outdir='/output';
    pplot = GEMplotter(['~/dogpack/apps/plasma/2d/twofluid/' models{1} outdir]);
    gplot = GEMplotter(['~/dogpack/apps/plasma/2d/twofluid/' models{2} outdir]);
    fnum = 1;
    pstate = read_state(pplot,fnum);
    gstate = read_state(gplot,fnum);
    pvar = sample_var_from_state(pplot,vars{i},pstate);
    gvar = sample_var_from_state(gplot,vars{i},gstate);
    plot(pplot,vars{i},160);
    pause
    dvar = pvar-gvar;
    %plot_scalar(pplot,dvar);
    display_frame(pplot,-1,dvar,['diff of ' vars{i}]);
    pause
  end
end

