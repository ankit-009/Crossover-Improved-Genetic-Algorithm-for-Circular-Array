function code = generateCode( obj )
%GENERATECODE Generate MATLAB code to recreate a GA options object
%   generateCode(gaoptimset ) is a string that can be evaluated to produce
%   the options structure passed in.

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:12 $

% read the list of properties from another file.
% this list is shared between several MATLAB files for consistency.
properties =  propertyList;

% make a default object. properties that match the values in the default
% will not be generated.
default = gaoptimset;

% first line
code = sprintf('options = gaoptimset;\n');

% for each property
for i = 1:length(properties)
    prop = properties{i};
    if(~isempty(prop)) % the property list has blank lines, ignore them
        value = obj.(prop);
        if(~isequal(value,default.(prop))) % don't generate code for defaults.
            code = [code sprintf('options.%s = %s;\n',prop,value2RHS(value))];
        end
    end
end
