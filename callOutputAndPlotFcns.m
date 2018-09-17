function stop = callOutputAndPlotFcns(options, ...
    state, fname, localSolution, localRunIndex, funcCount, bestX, bestFval)
%CALLOUTPUTANDPLOTFCNS Call output and plot functions
%
%   STOP = CALLOUTPUTANDPLOTFCNS(OPTIONS, STATE, FNAME, LOCALSOLUTION,
%   POINTINDEX, FUNCCOUNT, BESTX, BESTFVAL) calls the output and plot
%   functions specified in the solver OPTIONS structure. The optimValues
%   structure passed to each function will contain the fields
%   localSolution, pointIndex, funcCount, bestX and bestFval.
%
%   Note that this function assumes that the output and plot functions are
%   in options.OutputFcns and options.PlotFcns. Furthermore, these fields
%   are assumed to be either empty or cell arrays.
%
%   See also GADSPLOT

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:22 $

% Create optimValues structure. Lower case field names used to be
% consistent with the Simulated Annealing and Optimization toolbox solvers.
optimValues.localsolution = localSolution;
optimValues.localrunindex = localRunIndex;
optimValues.funccount = funcCount;
optimValues.bestx = bestX;
optimValues.bestfval = bestFval;

% Call all output functions
stop = false;
if ~isempty(options.OutputFcns)
    switch state
        case {'iter','init'}
            stop = i_callAllOutputFcns(options.OutputFcns, optimValues, state);
        case 'done'
            i_callAllOutputFcns(options.OutputFcns, optimValues, state);
        otherwise
            error(message('globaloptim:callOutputAndPlotFcns:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end

% Call plot functions
if ~isempty(options.PlotFcns)
    switch state
        case {'iter','init'}
            stop = gadsplot(options, optimValues, state, fname) || stop;
        case 'done'
            gadsplot(options, optimValues, state, fname);
        otherwise
            error(message('globaloptim:callOutputAndPlotFcns:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end

function stop = i_callAllOutputFcns(OutputFcn, optimValues, state)

% Call each output function
stop = false(length(OutputFcn),1);
for i = 1:length(OutputFcn)
    stop(i) = feval(OutputFcn{i},optimValues,state);
end

% If any stop(i) is true we set the stop to true
stop = any(stop);
