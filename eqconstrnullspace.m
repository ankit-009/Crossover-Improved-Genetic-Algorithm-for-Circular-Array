function NullSpace = eqconstrnullspace(A,nDims)
%EQCONSTRNULLSPACE is private to PATTERNSEARCH
% Function that returns the null space of the input matrix A .
% nDims is the dimensionality of the optimization space.
%
% Private to PATTERNSEARCH
%
% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:33 $
%

% One can also do this by finding the non trivial orthonormal basis from QR itself
% but it currently causes scaling issues in some of the test problems. 
% The NULL function uses SVD which is less efficient  and more appropriate if A 
% might not be full row rank, but rank(A) = size(A,1) is guaranteed by pre processing.

numEqConstr = size(A,1);
if ~isempty(A)
    % Find Null Space of equality constraints
    [Q,~] = qr(A',0);
    if (numEqConstr <= nDims)
		NullSpace = eye(nDims) - Q*Q';
    else
        error(message('globaloptim:eqconstrnullspace:Infeasible')); 
    end
else
    NullSpace = eye(nDims);
end

