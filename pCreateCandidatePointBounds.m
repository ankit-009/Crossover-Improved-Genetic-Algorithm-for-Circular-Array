function [candPtLB, candPtUB] = pCreateCandidatePointBounds(LB, UB)
%PCREATECANDIDATEPOINTBOUNDS Create candidate point bounds.
%
%   [CANDPTLB, CANDPTUB] = PCREATECANDIDATEPOINTBOUNDS(LB, UB) creates a
%   set of finite candidate point bounds to be used to generate a set of
%   candidate points for GLOBALSEARCHNLP.
%
%   Private to Global Optimization Toolbox

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:38 $

% First, shift the bounds by 1 unit to stop initial reference set for
% GLOBALSEARCHNLP containing the origin. We have found that this causes
% problems during testing.
ARTIFICIAL_LOWER_BOUND = -1e4 + 1;
ARTIFICIAL_UPPER_BOUND = 1e4 + 1;
candPtLB = LB;
candPtUB = UB;

% Deal with lower bounds. Replace non-finite elements with
% ARTIFICIAL_LOWER_BOUND. 
idxReplaceLB = ~isfinite(candPtLB);
candPtLB(idxReplaceLB) = ARTIFICIAL_LOWER_BOUND;

% Deal with lower bounds. Replace non-finite elements with
% ARTIFICIAL_UPPER_BOUND. 
idxReplaceUB = ~isfinite(candPtUB);
candPtUB(idxReplaceUB) = ARTIFICIAL_UPPER_BOUND;

% Deal with the case where the lower bound is now greater than the upper
% bound
idxLBgtUB = candPtLB > candPtUB;

% Case 1: Caused by replacement of lower bound
idxAddLB = idxReplaceLB & idxLBgtUB;
candPtLB(idxAddLB) = candPtUB(idxAddLB) + ARTIFICIAL_LOWER_BOUND;

% Case 2: Caused by replacement of upper bound
idxAddUB = idxReplaceUB & idxLBgtUB;
candPtUB(idxAddUB) = candPtLB(idxAddUB) + ARTIFICIAL_UPPER_BOUND;

% Note that LB cannot exceed UB if both bounds are replaced. Also we have
% already checked that LB <= UB for the finite bounds.

