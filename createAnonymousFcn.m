function fcn_handle = createAnonymousFcn(fcn,FcnArgs)
%CREATEFUNCTIONHANDLE create an anonymous function handle
%
% fcn: A function handle 
% args: A cell array of extra arguments to user's objective/constraint
% function

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $Date: 2012/08/21 00:23:29 $

fcn_handle = @(x) fcn(x,FcnArgs{:});
