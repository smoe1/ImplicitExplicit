
% could implement this with a call to read_config_file
%
function [gamma,mass_ratio,temp_ratio,cs_light,B_0,BCs,domain_scaling,cc,clean_E_field] ...
  = read_param(outputdir)

  readfile=[outputdir '/param.data'];

  % alternative implementation/invocation
  %
  global gamma mass_ratio temp_ratio cs_light B_0 BCs domain_scaling ...
     cc clean_E_field;
  domain_scaling=1;
  read_config_file(readfile);

 if(0) % commenting out
  fid=fopen(readfile,'r');
  if(fid==-1)
     error(['could not open file ' readfile]);
  end
  
  fgetl(fid); % model parameters
  fgetl(fid); %
  gamma         = get_num(fid, 'gamma' );
  mass_ratio    = get_num(fid, 'mass_ratio' );
  temp_ratio    = get_num(fid, 'temp_ratio' );
  cs_light      = get_num(fid, 'cs_light' );
  fgetl(fid); %
  fgetl(fid); % problem parameters (ICs and BCs)
  fgetl(fid); %
  B_0           = get_num(fid, 'B_0' );
  BCs           = get_num(fid, 'BCs' );
  domain_scaling= get_num(fid, 'domain_scaling' );
  fgetl(fid); %
  fgetl(fid); % method parameters
  fgetl(fid); %
  cc            = get_num(fid, 'cc' );
  clean_E_field = get_num(fid, 'clean_E_field' );
 end
end
