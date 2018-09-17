function options = initialPopulationCheck(options)
%initialPopulationCheck Validates initial population

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:17 $

% Make a *guess* at numberOfObjectives 
if ~isempty(options.InitialScores)
    numberOfObjectives = size(options.InitialScores,2);
else
    numberOfObjectives = 1; % safe for testing
end

% Scores for single objective should be a column vector
if numberOfObjectives == 1
   options.InitialScores = options.InitialScores(:);
end
size_initPop  = size(options.InitialPopulation,1);
[len_initScore,numObj] = size(options.InitialScores);

% No tests if initial pop and scores are empty
if size_initPop == 0 && len_initScore == 0
    return;
end

popSize  = sum(options.PopulationSize);

if size_initPop > popSize
    warning(message('globaloptim:initialPopulationCheck:initPopLength'));
    options.InitialPopulation(popSize+1:size_initPop,:) = [];
    size_initPop  = size(options.InitialPopulation);
end

if len_initScore > popSize
    warning(message('globaloptim:initialPopulationCheck:initScoreLength'));
    options.InitialScores(popSize+1:len_initScore,:) = [];
    len_initScore = size(options.InitialScores,1);
end

if len_initScore > size_initPop
    warning(message('globaloptim:initialPopulationCheck:initScoreAndPopLength'));
end

