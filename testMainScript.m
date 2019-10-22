%% Set the options
clear
clc
problemSize = 30; % Number of decision variables
nsga_opt = nsgaiiOptions(problemSize);

%% initial population
Population = generateInitialpopulation(nsga_opt);

%% evaluate population
for ii = 1:nsga_opt.PopulationSize
    f(ii,:) = testObjectiveFunctions(Population(ii,:), 'ZDT1');    
end

%% Non Dominating sorting
[rank, accum] = nonDominatedSorting(f, 'Goldberg');
%[rank, accum] = nonDominatedSorting(f, 'Fonseca');
plotParetoFronts(f, rank, 'Goldberg');
%% Selection
SelectedParents = parentSelection(rank, accum, nsga_opt);
%% Crossover 
xPopulation = crossoverParents(Population, SelectedParents, rank, accum, nsga_opt);
%% Mutation
xPopulation = mutation(xPopulation, nsga_opt);
%% create inline function
f = @(x)testObjectiveFunctions(x, 'ZDT3');
%% Run NSGAII
ParetoSolutions = nsgaii(f, nsga_opt);