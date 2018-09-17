function [handle,args] = functionHandleOrCell(property,value)
%functionHandleOrCell A function Handle or a cell array starting with a function
%handle.

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:38 $

[handle,args] = isFcn(value);

if ~isempty(handle)
    return
elseif strcmp(property,'NonconFcn')
    error(message('globaloptim:functionHandleOrCell:needFunctionHandleConstr'));
else
    error(message('globaloptim:functionHandleOrCell:needFunctionHandle', property));
end

