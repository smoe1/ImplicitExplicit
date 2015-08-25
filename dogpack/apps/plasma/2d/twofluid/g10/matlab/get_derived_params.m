
function [one_over_epsilon,elc_mass,ion_mass] = get_derived_params(cs_light,mass_ratio);
  one_over_epsilon = cs_light*cs_light;
  elc_mass = 1./(mass_ratio+1.);
  ion_mass = mass_ratio*elc_mass;
end
