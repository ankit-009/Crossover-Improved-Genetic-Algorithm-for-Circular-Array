function D = distancepoints(X, Y)
%DISTANCEPOINTS Return matrix of distances between points.
%
%  D = DISTANCEPOINTS(X, Y) generates a matrix containing the distance
%  between each point in X and Y.
%
%  D = DISTANCEPOINTS(X) generates a matrix containing the distance between
%  each point in X and every other point in X.
%
%  Notes:-
%  1) This function and its corresponding MEX function perform no error 
%  checking. 
%  
%  2) It is assumed that X and Y satisfy the following:
%     a. X and Y are both double arrays with real elements.
%     b. Each column of X and Y is a point.
%     c. Number of rows of X and Y are identical.
%
%  3) We are aware that other functions that calculate distance in MATLAB,
%  e.g. pdist in Statistics Toolbox and xregdistancepoints in MBC, assume
%  that each row of X and Y is a point. However, given that distancepoints
%  is a private function for GlobalSearch and MultiStart there isn't a
%  requirement on distancepoints to be consistent with these other distance
%  calculation functions. 

%   Copyright 2009-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:32 $

% Calculation performed in MEX function
if nargin==1
    D = mx_distancepoints(X, []);
    % In one input case, the matrix is symmetric with zeros on the
    % diagonal. The mex function just calculates the lower triangular
    % portion of D, we construct the rest here. 
    D = D + D';
else
    D = mx_distancepoints(X, Y);
end

