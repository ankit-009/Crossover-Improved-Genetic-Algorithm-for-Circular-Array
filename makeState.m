function state = makeState(GenomeLength,FitnessFcn,Iterate,type,options)
%MAKESTATE Create an initial population and fitness scores
%   state = makeState(GenomeLength,FitnessFcn,options) Creates an initial 
%   state structure using the information in the options structure.

%   Copyright 2003-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2012/08/21 00:24:29 $

% A variety of data used in various places
state.Generation = 0;		% current generation counter
state.StartTime = cputime;	% start time
state.StopFlag = []; 		% reason for termination
state.LastImprovement = 1;	% generation stall counter
state.LastImprovementTime = state.StartTime;	% time stall counter
state.Best = [];            % best score in every generation
state.how = '';
state.FunEval = 0;
popSize = sum(options.PopulationSize);
state.Expectation = zeros(popSize,1);  % expection of individuals
state.Selection = zeros(popSize,1);       % selection indices

% If InitialPopulation is partly empty we will use the creation function to
% generate population (CreationFcn can utilize InitialPopulation)
if sum(options.PopulationSize) ~= size(options.InitialPopulation,1)
    % If initial population is empty and we know one feasible individual then
    % assign it to InitialPopulation
    if ~isempty(Iterate.x) && isempty(options.InitialPopulation)
        options.InitialPopulation(1,:) = Iterate.x';
    end
      state.Population = feval(options.CreationFcn,GenomeLength,FitnessFcn,options,options.CreationFcnArgs{:});
else % the initial pop was passed in!
    state.Population = options.InitialPopulation;
end

% Calculate score for state.Population
totalPopulation = size(state.Population,1);
initScoreProvided = length(options.InitialScores);
if totalPopulation ~= initScoreProvided
   individualsToEvaluate = totalPopulation - initScoreProvided;
   state.Score  =  zeros(totalPopulation,1);
   
    if initScoreProvided > 0
        state.Score(1:initScoreProvided) = options.InitialScores(:);
    end
    % Score each member of the population
    if strcmpi(options.Vectorized, 'off')
        try
            firstMemberScore = FitnessFcn(state.Population(initScoreProvided+1,:));
        catch userFcn_ME
            gads_ME = MException('globaloptim:makestate:fitnessFcnFailed', ...
                'Failure in initial user-supplied fitness function evaluation. GA cannot continue.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end
        % User-provided fitness function should return a scalar
        if numel(firstMemberScore) ~= 1
            error(message('globaloptim:makestate:fitnessCheck'));
        end        
        Score = fcnvectorizer(state.Population(initScoreProvided+2:end,:),FitnessFcn,1,options.SerialUserFcn);
        % Concatenate the score of the first member of the population to the rest
        Score = [firstMemberScore; Score];
        state.Score(initScoreProvided+1:end) =  Score(:);
    else
        try
            Score = FitnessFcn(state.Population(initScoreProvided+1:end,:));
        catch userFcn_ME
            gads_ME = MException('globaloptim:makestate:fitnessFcnFailed', ...
                'Failure in initial user-supplied fitness function evaluation. GA cannot continue.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end
        if numel(Score) ~= individualsToEvaluate
            error(message('globaloptim:makestate:fitnessCheckVectorized', ...
                'Vectorized','on'));
        end
        state.Score(initScoreProvided+1:end) =  Score(:);
    end
    state.FunEval = individualsToEvaluate;          % number of function evaluations
else
    state.Score = options.InitialScores;
    state.FunEval = 0;          % number of function evaluations
end
% Make sure score is a column vector
state.Score = state.Score(:);

% Partial population is allowed for 'doubleVector' and 'bitString'
% population type so make population of appropriate size
if ~strcmpi(options.PopulationType,'custom')
    lens = size(state.Population,1);
    npop = sum(options.PopulationSize);
    if npop > lens
        population = zeros(npop,GenomeLength);
        population(1:lens,:) = state.Population;
        population(lens+1:end,:) = repmat(state.Population(end,:),(npop-lens),1);
        scores = zeros(npop,1);
        scores(1:lens) = state.Score;
        scores(lens+1:end) = repmat(state.Score(end),(npop-lens),1);
        state.Population = population;
        state.Score = scores;
    else
        state.Population(npop+1:end,:) = [];
        state.Score(npop+1:end) = [];
    end
end

