function innerOptions = psAugReset(options,optimState,psAugParam)
%PSAUGRESET Reset some variables before solving a new sub-problem
%   Private to PATTERNSEARCH.

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:55 $

% Reset values before inner iteration
innerOptions = options;
% Setting maxiter to Inf is okay because maxfuneval is still finite
innerOptions.MaxIter = Inf; %200*numberOfVariables;
innerOptions.TolX = psAugParam.currentTolMesh; 
innerOptions.TolFun = psAugParam.currentTolMesh;
innerOptions.TolMesh = psAugParam.currentTolMesh;
innerOptions.MaxFunEvals = options.MaxFunEvals - optimState.FunEval;
innerOptions.TimeLimit = options.TimeLimit - (cputime-optimState.StartTime);
innerOptions.Verbosity = 0;
% 'Dynamic' options for ScaleMesh means that scale is updated in every
% iteration (undocumented and used internally by nonlinear constrained solver)
if strcmp(options.ScaleMesh,'on')
   innerOptions.ScaleMesh = 'dynamic';
end

