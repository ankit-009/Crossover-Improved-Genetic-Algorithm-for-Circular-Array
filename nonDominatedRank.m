function nondominatedRank = nonDominatedRank(score,nParent)
%nonDominatedRank Assigns rank to individuals in 'score'. 
%   The optional argument 'nParent' will limit rank assignment to only
%   'nParent' individuals.

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:32 $


% Population size
popSize = size(score,1);
if nargin < 2
    nParent = popSize;
end
% Boolean to track if individuals are ranked
rankedIndiv = false(popSize,1);
% Initialization: rank of individuals (denotes which front)
nondominatedRank = inf(popSize,1);

numObj = size(score,2);
dominationMatrix = false(popSize);
% First test for domination is to check which points have better function
% values in at least one of the objectives
for count = 1:numObj
    dominationMatrix = dominationMatrix | bsxfun(@lt,score(:,count),score(:,count)');
end
% Now, check to see if those points that pass the first test, if they are
% at least as good as others in all the objectives
for count = 1:numObj
    dominationMatrix = dominationMatrix & bsxfun(@le,score(:,count),score(:,count)');
end
% We will do the test along the column
dominationMatrix = ~dominationMatrix;

% At this point, we have the domination matrix that may look like this
% (example only). 
    %     p1  p2  p3  p4  p5
    % p1  1   1   0   1   1
    %
    % p2  1   1   0   1   1
    %
    % p3  0   1   1   1   1
    %
    % p4  0   0   1   1   0
    %
    % p5  1   1   0   1   1

% In this matrix, if (i,j) entry is '1' then jth individual dominates ith
% individual else it does not dominate ith individual. Also, if all the
% entries in a column are '1' then it is a non-dominated individual. In
% this example we have 5 points. From this matrix we can tell that 4th
% point is non-dominated since 4th column has all '1'. In the first
% iteration (rank = 1) p4 will be chosen as non-dominated individual.

rankCounter = 1;
while ~all(rankedIndiv) && nnz(isfinite(nondominatedRank)) <= nParent
    dominates = all(dominationMatrix);
    nondominatedRank(dominates) = rankCounter;
    rankCounter = rankCounter + 1;
%=========================================================================
    % We want to remove p4 from the population but we don't want to change
    % the matrix size. Instead, we turn all the entris in the 4th row to
    % '1' which means that p4 is dominated by every other individuals
    % (effectively).
    dominationMatrix(dominates,:) = true;
        
    %     p1  p2  p3  p4  p5
    % p1  0   1   0   1   1
    %
    % p2  1   1   0   1   1
    %
    % p3  0   1   1   1   1
    %
    % p4  1   1   1   1   1   
    %
    % p5  1   1   0   1   1
    
    % By this action, p2 and p5 are now non-dominated individuals
    % (excluding p3) in the remaining pool.
    
%==========================================================================
    % Next, make sure we don't pick the same individuals again. For this to
    % happen, we turn the (4,4) element to '0'.
    dominationMatrix(dominates,dominates) = false;
    
    %     p1  p2  p3  p4  p5
    % p1  1   1   0   1   1
    %
    % p2  1   1   0   1   1
    %
    % p3  0   1   1   1   1
    %
    % p4  1   1   1   0   1   
    %
    % p5  1   1   0   1   1
    
    % In the next iteration (rank = 2) p5 and p2 will be chosen but not p4. 
%==========================================================================
    
    rankedIndiv(dominates) = true;
end


% For the above example, it will take 4 iterations to find exactly 4 pareto
% fronts (with rank = 1, 2, 3, and 4). At the end of the 2nd iteration, the
% matrix will look like this (two fronts found)
    %     p1  p2  p3  p4  p5
    % p1  1   1   1   1   1
    %
    % p2  1   0   1   1   1
    %
    % p3  0   1   1   1   1
    %
    % p4  1   1   1   0   1   
    %
    % p5  1   1   1   1   0

% and, after the 3rd iteration it will look like this (3rd front will be
% found)

    %     p1  p2  p3  p4  p5
    % p1  1   1   1   1   1
    %
    % p2  1   0   1   1   1
    %
    % p3  1   1   0   1   1
    %
    % p4  1   1   1   0   1   
    %
    % p5  1   1   1   1   0

% and after the 4th iteration (last) the matrix will be not used anymore
% but it will look like this:

    %     p1  p2  p3  p4  p5
    % p1  0   1   1   1   1
    %
    % p2  1   0   1   1   1
    %
    % p3  1   1   0   1   1
    %
    % p4  1   1   1   0   1   
    %
    % p5  1   1   1   1   0
    
    