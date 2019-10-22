function ParetoSolutions = nsgaii(fun, options)

Population = generateInitialpopulation(options);
ParetoSolutions.x = [];
ParetoSolutions.f = [];

while options.Generation < options.MaxGeneration
    
    % evaluate population
    clear f
    for ii = 1:options.PopulationSize
        f(ii,:) = fun(Population(ii,:));
    end
    
    %% Non Dominating sorting
    f = [f;ParetoSolutions.f];
    Population = [Population; ParetoSolutions.x];
    [f, Population] = removeDuplicateSolutions(f, Population, options);
    
    [rank, accum] = nonDominatedSorting(f, options.DominatingSortingType);
    
    clf
    plotParetoFronts(f, rank, options);
    title(['Generation: ' num2str(options.Generation)]);
    drawnow
    
    ParetoSolutions.f = f(rank == 1,:);
    ParetoSolutions.x = Population(rank == 1,:);
    
    selectedIds = selectSurvivors(rank, accum, options);
    Population = Population(selectedIds,:);
    rank = rank(selectedIds,:);
    accum = accum(selectedIds,:);
    
    SelectedParents = parentSelection(rank, accum, options);
    
    Population = crossoverParents(Population, SelectedParents, rank, accum, options);
    
    Population = mutation(Population, options);
    
    options.Generation = options.Generation + 1;
end