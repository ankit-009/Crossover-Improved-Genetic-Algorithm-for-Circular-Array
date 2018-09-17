function range = checkPopulationInitRange(lowerBound,upperBound,range)
%checkPopulationInitRange performs check on bounds and PopInitRange and
%   make sure that they are consistent 

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:23:23 $

% Range is a 2 X nvars matrix
lowerRange = range(1,:)';
upperRange = range(2,:)';

meanLowerRange = mean(lowerRange(isfinite(lowerRange)));
meanUpperRange = mean(upperRange(isfinite(upperRange)));


lowerInf = isinf(lowerRange) & isinf(lowerBound);
upperInf = isinf(upperRange) & isinf(upperBound);

lowerRange(lowerInf) = meanLowerRange;
upperRange(upperInf) = meanUpperRange;

% When lower bound is greater than lower init range, choose lower bound
index = false(length(lowerRange),1);
index(~lowerInf) = lowerBound(~lowerInf) >= lowerRange(~lowerInf);
lowerRange(index) = lowerBound(index);

% When upper bound is greater than upper init range, choose upper init range
index = false(length(upperRange),1);
index(~upperInf) = upperBound(~upperInf) < upperRange(~upperInf);
upperRange(index) = upperBound(index);

% Make sure that range is consistent (lowerRange < upperRange)
% if not consistent we take the value from bounds if they are finite or from range
index  = (upperRange <= lowerRange) & isfinite(upperBound) ;
index1 = (upperRange <= lowerRange) & isinf(upperBound) ;
upperRange(index) = upperBound(index);
upperRange(index1) = lowerRange(index1) + 1; % Default range is [0 1]

index  = (upperRange <= lowerRange) & isfinite(lowerBound) ;
index1 = (upperRange <= lowerRange) & isinf(lowerBound) ;
lowerRange(index) = lowerBound(index);
lowerRange(index1) = upperRange(index1) - 1; % Default range is [0 1]

range = [lowerRange';upperRange'];
