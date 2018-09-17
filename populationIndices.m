function indices = populationIndices(population)
%POPULATIONINDICES Find the indices of each sub-population
%   indices = populationIndices(population); returns a 2 by n array
%   containing the locations of each subpopulation in the population array.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:45 $

% Renamed old file 'populationIndicies.m' to 'populationIndices.m' (note
% the typo in 'indices')

lengths = cumsum(population);
lengths = lengths(1:(end-1));
starts = [1, 1 + lengths];
ends = starts + population - 1;

indices = [starts;ends];