function positiveIntegerArray(property,value)
%positiveIntegerArray positive integer array

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:48 $

allValid = true;
for i = 1:numel(value)
    valid =  isreal(value(i)) && value(i) == floor(value(i)) && value(i) > 0;
    allValid = allValid && valid;
end
if(~valid)
    error(message('globaloptim:positiveIntegerArray:notPosIntegerArray', property));
end
