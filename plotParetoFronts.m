function plotParetoFronts(f, rank, options)

nRanks = max(rank);

hold on
switch options.DominatingSortingType
    case 'Goldberg'
        for ii = 1:nRanks
            fRank = f(rank == ii,:);
            [x, d] = sort(fRank(:,1));
            plot(x, fRank(d,2), '.:')
            if ~options.PlotAllFronts
                break;
            end
            text(x(1), fRank(d(1),2), num2str(ii), 'fontsize', 12);
        end
    case 'Fonseca'
        plot(f(:,1), f(:,2),'.')
        for ii = 1:length(rank)
            text(f(ii,1), f(ii,2), num2str(rank(ii)));
        end
        
end