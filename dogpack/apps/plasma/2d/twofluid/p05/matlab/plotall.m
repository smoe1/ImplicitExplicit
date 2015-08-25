
function plotall(outputdir_in)
  global outputdir;
  if(nargin>=1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end

  plot_recon(outputdir);
  %plotr(outputdir);
  plotout(outputdir);
end
