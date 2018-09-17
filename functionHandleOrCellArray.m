function [FUN,ARGS] = functionHandleOrCellArray(property,value)
%functionHandleOrCellArray A function Handle or a cell array starting with a
%   function handle.

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:23:39 $

FUN  = [];
ARGS = {};

% If a scalar  ~cell  is passed convert to cell (for clarity, not speed)
if ~iscell(value) && length(value) == 1
     value = {value};
end
fun_counter = 1;
% If value is an array of functions, it must be a cell array
for i = 1:numel(value)
    candidate = value(i);
    % If any element is also a cell array
    if iscell(candidate)
        if isempty(candidate{1})
            continue;
        end
        %Sometimes the variable 'candidate' might have nested cell array 
        %e.g. {{@outputfcn, p1,p2}} instead of just
        %{@outputfcn,p1,p2}. The following code gets rid of extra braces,
        %which are typically introduced by GUI import/export options.
        temp = candidate{1};
        while iscell(temp) && isscalar(temp)
            candidate = temp(1);
            temp = candidate{1};
        end
        [handle,args] = isFcn(candidate{:});
    else
        [handle,args] = isFcn(candidate);
    end
    if(~isempty(handle)) && (isa(handle,'inline') || isa(handle,'function_handle'))
        FUN{fun_counter} = handle;
        ARGS{fun_counter} = args;
        fun_counter = fun_counter + 1;
    else
        error(message('globaloptim:functionhandleorcellarray:needFunctionHandle', property));
    end
end
