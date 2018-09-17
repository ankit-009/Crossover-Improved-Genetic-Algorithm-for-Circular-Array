function options  = checkoptions(options,defaultopt,numberOfVariables)
%CHECKOPTIONS validates all PATTERNSEARCH options before they are used by
%   solver
%
%   private to pfminlcon, pfminbnd, and pfminunc.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:23:27 $

% Sanity check for the options structure
options = psoptimset(options);

options.Display = psoptimget(options,'Display',defaultopt,'fast');

% Define verbosity here (Later we can use options structure)
switch options.Display  
    case {'off','none'}
        options.Verbosity = 0;
    case 'final'
        options.Verbosity = 1;    
    case 'iter'
        options.Verbosity = 2;
    case 'diagnose'
        options.Verbosity = 3;
    otherwise
        options.Verbosity = 1;
end

% Retrieve options using PSOPTIMGET
options.MeshExpansion     = psoptimget(options,'MeshExpansion',defaultopt,'fast'); 
options.MeshContraction   = psoptimget(options,'MeshContraction',defaultopt,'fast'); 
options.CompleteSearch    = psoptimget(options,'CompleteSearch',defaultopt,'fast');
options.MeshAccelerator   = psoptimget(options,'MeshAccelerator',defaultopt,'fast');
options.TolMesh           = psoptimget(options,'TolMesh',defaultopt,'fast');
options.TolCon            = psoptimget(options,'TolCon',defaultopt,'fast');
options.MaxMeshSize       = psoptimget(options,'MaxMeshSize',defaultopt,'fast');
options.MaxIter           = psoptimget(options,'MaxIter',defaultopt,'fast');
options.MaxFunEvals       = psoptimget(options,'MaxFunEvals',defaultopt,'fast');
options.TimeLimit         = psoptimget(options,'TimeLimit',defaultopt,'fast');
options.TolBind           = psoptimget(options,'TolBind',defaultopt,'fast');
options.TolFun            = psoptimget(options,'TolFun',defaultopt,'fast');
options.TolX              = psoptimget(options,'TolX',defaultopt,'fast');
options.InitialMeshSize   = psoptimget(options,'InitialMeshSize',defaultopt,'fast');
options.PollMethod        = psoptimget(options,'PollMethod',defaultopt,'fast');
options.PollingOrder      = psoptimget(options,'PollingOrder',defaultopt,'fast');
options.CompletePoll      = psoptimget(options,'CompletePoll',defaultopt,'fast');
options.PlotInterval      = psoptimget(options,'PlotInterval',defaultopt,'fast');
options.Vectorized        = psoptimget(options,'Vectorized',defaultopt,'fast');
options.Cache             = psoptimget(options,'Cache',defaultopt,'fast');
options.CacheTol          = psoptimget(options,'CacheTol',defaultopt,'fast');
options.CacheSize         = psoptimget(options,'CacheSize',defaultopt,'fast');
options.ScaleMesh         = psoptimget(options,'ScaleMesh',defaultopt,'fast');
options.MeshRotate        = psoptimget(options,'MeshRotate',defaultopt,'fast');
options.InitialPenalty    = psoptimget(options,'InitialPenalty',defaultopt,'fast');
options.PenaltyFactor     = psoptimget(options,'PenaltyFactor',defaultopt,'fast');
options.UseParallel       = psoptimget(options,'UseParallel',defaultopt,'fast');
% These options will be stuffed in the structure later (after some
% processing)
outputFcns        = psoptimget(options,'OutputFcns',defaultopt,'fast');
plotFcns          = psoptimget(options,'PlotFcns',defaultopt,'fast');
searchFcn        = psoptimget(options,'SearchMethod',defaultopt,'fast');
           
% Modify some fields if they are not yet assigned
if ischar(options.MaxFunEvals) && isequal(lower(options.MaxFunEvals),'2000*numberofvariables')
        options.MaxFunEvals = 2000*numberOfVariables;
end
if ischar(options.MaxIter) && isequal(lower(options.MaxIter),'100*numberofvariables')
        options.MaxIter = 100*numberOfVariables;
end

options.MaxFunEvals  = floor(options.MaxFunEvals);
options.MaxIter = floor(options.MaxIter);

% If searchFcn is a cell array with additional arguments, handle them
if iscell(searchFcn)
    searchFcnArg = searchFcn(2:end);
    searchFcn = searchFcn{1};
else
    searchFcnArg = {};
end
% Search technique could be [], char, or function_handle
if isempty(searchFcn)
    searchFcn = [];
elseif isa(searchFcn,'function_handle')
    searchFcn = fcnchk(searchFcn);
    searchFcnString = func2str(searchFcn);
elseif ischar(searchFcn) 
    searchFcn = fcnchk(searchFcn);
    searchFcnString = func2str(searchFcn);
else
    error(message('globaloptim:checkoptions:invalidSearchMethod'));    
end

% Make sure that search method is different from poll method and not 'none'
if ~isempty(searchFcn) && any(strcmpi(searchFcnString,{options.PollMethod,'none'}))
   searchFcn = [];
end

% Only some choices can be strings (special case)
if isa(searchFcn,'function_handle') && any(strcmpi(searchFcnString,{'positivebasisnp1', 'positivebasis2n', ...
            'gpspositivebasisnp1','gpspositivebasis2n','madspositivebasisnp1', 'madspositivebasis2n', ...
            'gsspositivebasisnp1','gsspositivebasis2n'}))
        searchFcn = searchFcnString; % Convert to a string because these are not functions
end


options.SearchMethod = searchFcn;
options.SearchMethodArg = searchFcnArg;

% If options.MaxMeshSize is less than options.Meshsize (This should not happen)
if options.MaxMeshSize < options.InitialMeshSize
    warning(message('globaloptim:checkoptions:maxMeshSize'));
    options.InitialMeshSize = options.MaxMeshSize;
end

% It is NOT vectorized in these conditions
options.NotVectorizedPoll   = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompletePoll,'off')));
options.NotVectorizedSearch = (strcmpi(options.Vectorized,'off') || ...
    (strcmpi(options.Vectorized, 'on') && strcmpi(options.CompleteSearch,'off')));

% If using 2N basis or MADS RotatePattern has no effect.
if any(strcmpi(options.PollMethod,{'positivebasis2n','gpspositivebasis2n', ...
    'madspositivebasisnp1','madspositivebasis2n','gsspositivebasis2n'}))
    options.MeshRotate = 'off';
end
% Error checking on InitialPenalty and PenaltyFactor
if options.InitialPenalty < 1
    warning(message('globaloptim:checkoptions:smallPenalty'));
    options.InitialPenalty = defaultopt.InitialPenalty;
end
% Penalty factor to increase penalty
if options.PenaltyFactor <= 1
    warning(message('globaloptim:checkoptions:smallPenaltyFactor'));
    options.PenaltyFactor = defaultopt.PenaltyFactor;
end
% If outputFcns is a cell array with additional arguments, handle them
[options.OutputFcns,options.OutputFcnsArg] = functionHandleOrCellArray('OutputFcns',outputFcns);

 if isempty(options.OutputFcns)
     options.OutputTrue = false;
 else
     options.OutputTrue = true;
 end
 % If plotFcns is a cell array with additional arguments, handle them
[options.PlotFcns,options.PlotFcnsArgs] = functionHandleOrCellArray('PlotFcns',plotFcns);

 if isempty(options.PlotFcns)
     options.PlotTrue = false;
 else
     options.PlotTrue = true;
 end
 
% Test for valid strings
if ~isempty(options.UseParallel)
    stringSet('UseParallel',options.UseParallel,{'never','always'});
    options.SerialUserFcn = strcmpi(options.UseParallel,'never');
else
    options.SerialUserFcn = true;
end

