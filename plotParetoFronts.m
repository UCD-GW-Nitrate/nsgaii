function plotParetoFronts(f, rank, method)

nRanks = max(rank);

hold on
switch method
    case 'Goldberg'
        for ii = 1:nRanks
            fRank = f(rank == ii,:);
            [x, d] = sort(fRank(:,1));
            plot(x, fRank(d,2), '.:')
            text(x(1), fRank(d(1),2), num2str(ii), 'fontsize', 12);
        end
    case 'Fonseca'
        plot(f(:,1), f(:,2),'.')
        for ii = 1:length(rank)
            text(f(ii,1), f(ii,2), num2str(rank(ii)));
        end
        
end