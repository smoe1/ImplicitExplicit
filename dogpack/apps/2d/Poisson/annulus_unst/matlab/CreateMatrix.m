function CreateMatrix(morder,showplots)
%CreateMatrix(morder, showplots)
%
% This is a single function that creates the matrix needed for the Poisson
% solve.
%
% Input parameters:
%
%    morder - order of the method.  This needs to match two items: subfactor in
%    the mesh generator as well as space_order in parameters.ini.
%
%    showplots = 1:   Show the plots (defualt value)
%              = 0:   Do not show the plots.
%

addpath('../../lib/matlab');
if( nargin < 2 )
    showplots=1;
end
bctype=1;
CreateMatrixGlobal(morder,bctype,showplots);
