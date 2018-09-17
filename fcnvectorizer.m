function y = fcnvectorizer(pop,fun,numObj,SerialUserFcn)
%FCNVECTORIZER is a utility function used for scalar fitness functions.

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:35 $

try
    popSize = size(pop,1);
    y = zeros(popSize,numObj);
    % Use for if SerialUserFcn is 'true'; the if-and-else part should
    % be exactly same other than parfor-for syntax difference.
    if SerialUserFcn
        for i = 1:popSize
            y(i,:) = feval(fun,(pop(i,:)));
        end
    else
        parfor (i = 1:popSize)
            y(i,:) = feval(fun,(pop(i,:)));
        end
    end
catch userFcn_ME
    if iscell(pop) % if it is a custom PopulationType
        y = fun(pop);
    else
         gads_ME = MException('globaloptim:fcnvectorizer:fitnessEvaluation', ...
             'Failure in user-supplied fitness function evaluation. GA cannot continue.');
         userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
end
