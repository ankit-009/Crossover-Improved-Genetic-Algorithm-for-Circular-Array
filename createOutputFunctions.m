function OPTIONS = createOutputFunctions(OPTIONS, timeStore)
%CREATEOUTPUTFUNCTIONS Create output functions for output structure storage
%and MaxTime limit.
%   CREATEOUTPUTFUNCTIONS creates two output functions and appends them to
%   the cell array of output functions provided in OPTIONS. First of these
%   output functions is a container to store output information. This will
%   allow us to retrieve part of the output structure if an error occurs
%   when running the local solver. The second output function is to allow a
%   user to set a maximum clock time for the global algorithm and will be
%   created only if timeStore.MaxTime is finite. OPTIONS at the output will
%   contain two (one if timeStore.MaxTime is infinite) output functions
%   under the 'OutputFcn' field that are appended to the existing list of
%   output functions. timeStore is a structure with fields 'maxTime' and
%   'startTime'. 
%
%   See also FMULTISTART, GLOBALSEARCHNLP

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:23:30 $

% Create output store
outStoreCont = outputStoreContainer;

% Create output function to update output information
hOutStore = @(x, optimValues, state)outstorefcn(x, optimValues, state, outStoreCont);

% Get existing output functions and convert to a cell array
outFcns = optimget(OPTIONS, 'OutputFcn');
if isempty(outFcns)
    outFcns = cell(0, 0);
end
if ~iscell(outFcns)
    outFcns = {outFcns};
end

% Output functions are now a cell, but can be any dimension. Make them a
% column vector.
outFcns = outFcns(:);

% Create output function to allow a user to set a maximum time for the
% algorithm, if maximum time is finite. Concatenate this function with the
% output store function.
if isinf(timeStore.maxTime)
    hGSOutFcns = {hOutStore};
else
    hTimeStopper = @(x, optimValues, state) timestorefcn(x, optimValues, ...
        state, timeStore);
    hGSOutFcns = {hTimeStopper; hOutStore};
end

% Append output functions and set in OPTIONS. 
outFcns = [outFcns; hGSOutFcns];    
OPTIONS = optimset(OPTIONS, 'OutputFcn', outFcns);
