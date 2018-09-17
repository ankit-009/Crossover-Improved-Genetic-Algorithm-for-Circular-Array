function [optimState,nextIterate,MeshSize,EXITFLAG,run] = ...
    getinitial(Iterate,numberOfVariables,neqcstr,lb,ub,options)
%GETINITIAL is private to pfminlcon, pfminbnd and pfminunc.

%   Copyright 2003-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:14 $

% Initialization
optimState.Iter = 0;
optimState.FunEval = 1;
optimState.infMessage  = '';
optimState.stopPlot = false;
optimState.stopOutput = false;
optimState.deltaF = NaN;
optimState.deltaX = NaN;
optimState.Successdir = [];
optimState.how = ' ';
optimState.MeshCont = options.MeshContraction;
optimState.scale = ones(numberOfVariables,1);
EXITFLAG = -1;
run = true;
MeshSize = options.InitialMeshSize;

% Calculate scale
if any(strcmpi(options.ScaleMesh,{'dynamic','on'})) && ~neqcstr
    meanX = mean([Iterate.x],2);
    optimState.scale = logscale(lb,ub,meanX);
end

nextIterate = Iterate;
optimState.StartTime = cputime;
