
% could implement this with a call to read_config_file
%
function [gamma,mass_ratio,temp_ratio,cs_light,B_0,BCs,cc,clean_E_field] ...
  = read_param(outputdir)

  readfile=[outputdir '/param.data'];

  global gamma mass_ratio temp_ratio cs_light B_0 BCs domain_scaling  ...
     cc clean_E_field enforced_symmetry;
  global gamma B_0 BCs domain_scaling cc enforced_symmetry;
  % set default values
  domain_scaling=1;
  enforced_symmetry=3;
  % read the values
  read_config_file(readfile);

end
