function state = gamultiobjMakeState(GenomeLength,FitnessFcn,type,options)
%gamultiobjMakeState Create an initial population and fitness scores

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:08 $

% A variety of data used in various places
state.Generation = 0;		% current generation counter
state.StartTime = cputime;	% start time
state.StopFlag = []; 		% reason for termination
state.FunEval = 0;
state.Selection = [];       % selection indices
popSize = sum(options.PopulationSize);

% If InitialPopulation is partly empty we will use the creation function to
% generate population (CreationFcn can utilize InitialPopulation)
if sum(options.PopulationSize) ~= size(options.InitialPopulation,1)
      state.Population = feval(options.CreationFcn,GenomeLength,FitnessFcn,options,options.CreationFcnArgs{:});
else % Initial population was passed in
    state.Population = options.InitialPopulation;
end

% Evaluate fitness function to get the number of objectives
try
    Score = FitnessFcn(state.Population(1,:));
catch userFcn_ME
    gads_ME = MException('globaloptim:gamultiobjMakeState:fitnessCheck', ...
        'Failure in initial user-supplied fitness function evaluation. GAMULTIOBJ cannot continue.');
    userFcn_ME = addCause(userFcn_ME,gads_ME);
    rethrow(userFcn_ME)
end
state.FunEval = state.FunEval + 1;
Score = Score(:)';
numObj = numel(Score);
% Size of InitialScore and Score should match
if ~isempty(options.InitialScores) && numObj ~= size(options.InitialScores,2)
    error(message('globaloptim:gamultiobjMakeState:initScoreSize','size(InitialScore,2)'));
end

% Calculate score for state.Population
totalPopulation = size(state.Population,1);
initScoreProvided = size(options.InitialScores,1);
if totalPopulation ~= initScoreProvided
    individualsToEvaluate = totalPopulation - initScoreProvided;
    state.Score  =  zeros(totalPopulation,numObj);

    if initScoreProvided > 0
        state.Score(1:initScoreProvided,:) = options.InitialScores(1:initScoreProvided,:);
    end
    % Score remaining members of the population
    if strcmpi(options.Vectorized, 'off')
        Score = fcnvectorizer(state.Population(initScoreProvided+1:end,:),FitnessFcn,numObj,options.SerialUserFcn);
        state.Score(initScoreProvided+1:end,:) =  Score;
    else
        Score = FitnessFcn(state.Population(initScoreProvided+1:end,:));
        if size(Score,1) ~= individualsToEvaluate
           error(message('globaloptim:gamultiobjMakeState:fitnessVectorizedCheck', ... 
                'Vectorized','on'));
        end
        state.Score(initScoreProvided+1:end,:) =  Score;
    end
    state.FunEval = state.FunEval+ individualsToEvaluate;          % number of function evaluations
else
    state.Score = options.InitialScores;
end

% Partial population is allowed for 'doubleVector' and 'bitString'
% population type so make population of appropriate size
if ~strcmpi(options.PopulationType,'custom')
    lens = size(state.Population,1);
    npop = sum(options.PopulationSize);
    if npop > lens
        population = zeros(npop,GenomeLength);
        population(1:lens,:) = state.Population;
        population(lens+1:end,:) = repmat(state.Population(end,:),(npop-lens),1);
        scores = zeros(npop,numObj);
        scores(1:lens,:) = state.Score;
        scores(lens+1:end,:) = repmat(state.Score(end,:),(npop-lens),1);
        state.Population = population;
        state.Score = scores;
    else
        state.Population(npop+1:end,:) = [];
        state.Score(npop+1:end,:) = [];
    end
end

% Get the rank and Distance measure of the population
[state.Population,state.Score,state.Rank,state.Distance] = ...
    rankAndDistance(state.Population,state.Score,options,popSize);

% One step GA (selection, crossover and then mutation)
% How many crossover offspring will there be from each source?
nXoverKids = round(options.CrossoverFraction * popSize);
nMutateKids = popSize - nXoverKids;
% how many parents will we need to complete the population?
nParents = 2 * nXoverKids + nMutateKids;
% Selection.
parents = feval(options.SelectionFcn,[state.Rank,state.Distance],nParents,options,options.SelectionFcnArgs{:});
% If selection function does not return the correct parents, error here
if (length(parents) ~= nParents) || min(parents) <= 0 || max(parents) > popSize
   error(message('globaloptim:gamultiobjMakeState:selectionFcnError', nParents, popSize));
end
% Shuffle to prevent locality effects. 
parents = parents(randperm(length(parents)));

xoverKids  = feval(options.CrossoverFcn,parents(1:(2 * nXoverKids)),options,GenomeLength, ...
     FitnessFcn,state.Score,state.Population,options.CrossoverFcnArgs{:});
mutateKids = feval(options.MutationFcn,parents((1 + 2 * nXoverKids):end),options,GenomeLength, ...
    FitnessFcn,state,state.Score,state.Population,options.MutationFcnArgs{:});

nextPopulation = [xoverKids ; mutateKids ];

if strcmpi(options.Vectorized, 'off') 
    nextScore = fcnvectorizer(nextPopulation,FitnessFcn,numObj,options.SerialUserFcn);
else
    nextScore = FitnessFcn(nextPopulation);
end
 state.FunEval =  state.FunEval + size(nextPopulation,1);
 
% Repeat for each sub-population (element of the PopulationSize vector)
offset = 0;
totalPop = options.PopulationSize;
for pop = 1:length(totalPop)
    populationSize =  totalPop(pop);
    subpopIndex = 1 + (offset:(offset + populationSize - 1));

    % Combine old and new population
    population = [state.Population(subpopIndex,:);nextPopulation(subpopIndex,:)];
    score_old = state.Score(subpopIndex,:);
    score = [state.Score(subpopIndex,:);nextScore(subpopIndex,:)];
    popSize = numel(subpopIndex);

    % Calculate the rank and distance measure of each individual 
    [state.Population(subpopIndex,:),state.Score(subpopIndex,:), ...
        state.Rank(subpopIndex,:),state.Distance(subpopIndex,:)] = ...
        rankAndDistance(population,score,options,popSize);
    offset = offset + populationSize;

    % Calculate average distance and spread
    [state.AverageDistance(pop), state.Spread(state.Generation+1,pop)] = ...
    distanceAndSpread(state.Distance(subpopIndex,:),state.Rank(subpopIndex,:), ...
    state.Score(subpopIndex,:),score_old);
end

