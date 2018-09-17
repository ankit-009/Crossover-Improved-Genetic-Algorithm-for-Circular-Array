function nonNegScalar(property,value)
%nonNegScalar any scalar >= 0    

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:34 $

valid =  isreal(value) && isscalar(value) && (value >= 0);
if(~valid)
    error(message('globaloptim:nonNegScalar:notNonNegScalar', property));
end
