function indices = populationIndicies(population)
%POPULATIONINDICIES Find the indices of each sub-population
%   POPULATIONINDICIES has been replaced by POPULATIONINDICES, and will be 
%   removed in a later version
%   indices = populationIndicies(population); returns a 2 by n array
%   containing the locations of each subpopulation in the population array.
%
%   Private to GA

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:46 $


warning(message('globaloptim:populationIndicies:obsoleteFcn'));
indices = populationIndices(population);

