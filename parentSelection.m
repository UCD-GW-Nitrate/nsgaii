function SelectedParents = parentSelection(rank, accum, options)

switch options.Crossover.Method
    case {'Uniform', 'Intermediate'}
        SelectedParents = nan(options.PopulationSize,1);
    case 'Heuristic'
        SelectedParents = nan(2*options.PopulationSize,1);
end

switch options.Selection.Method
    case 'Tournament'
        for ii = 1:length(SelectedParents)
            while true
                tourid = randperm(options.PopulationSize ,options.Selection.TourSize);
                minRank = min(rank(tourid));
                idsWithMinRank = find(rank(tourid) == minRank);
                if length(idsWithMinRank) == 1
                    SelectedParents(ii,1) = tourid(idsWithMinRank);
                else
                    [~, d] = sort(accum(tourid(idsWithMinRank)));
                    SelectedParents(ii,1) = tourid(idsWithMinRank(d(end)));
                end
                if rem(ii,2) == 0
                    if SelectedParents(ii,1) ~= SelectedParents(ii-1,1)
                        break;
                    end
                else
                    break;
                end
            end
        end
end