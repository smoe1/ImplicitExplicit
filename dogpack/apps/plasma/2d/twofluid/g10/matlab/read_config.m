% read the parameters from a file in the format
%
%  1.2    variable_name[1] = comment 1
%  2      variable_name[2] % comment 2
%
%  str    variable_name_3  =comment 3
%
global outputdir;
if(isempty(outputdir))
    outputdir='output';
end

% define derived parameters
%
global mass_ratio;
global elc_mass ion_mass;

appConfig = [outputdir '/param.data'];
read_config_file(appConfig);
libConfig = [outputdir '/dogpack.data'];
read_config_file(libConfig);

elc_mass = 1./(mass_ratio+1.);
ion_mass = mass_ratio*elc_mass;
