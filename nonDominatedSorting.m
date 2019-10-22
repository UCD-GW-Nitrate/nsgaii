function [rank, accum] = nonDominatedSorting(f, method)

populationSize = size(f,1);
rank = nan(populationSize,1);
switch method
    case 'Goldberg'
        poolID = [1:populationSize]';
        currentRank = 1;
        while ~isempty(poolID)
            dltID = [];
            for ii = 1:length(poolID)
                test = f(poolID,1) < f(poolID(ii),1);
                for jj = 2:size(f,2)
                    test = test & f(poolID,jj) < f(poolID(ii),jj);
                end
                nDominatingSolutions = sum(test);
                if nDominatingSolutions == 0
                    rank(poolID(ii),1) = currentRank;
                    dltID = [dltID;ii];
                end
            end
            poolID(dltID,:) = [];
            currentRank = currentRank + 1;
        end
    case 'Fonseca'
        for ii = 1:populationSize
            test = f(:,1) < f(ii,1);
            for jj = 1:size(f,2)
                test = test & f(:,jj) < f(ii,jj);
            end
            rank(ii,1) = sum(test);
        end
end

nRanks = max(rank);
accum = zeros(populationSize,1);
fmin = min(f,[], 1);
fmax = max(f,[], 1);
for ii = 1:nRanks
    id = find(rank == ii);
    for kk = 1:size(f,2)
        [c, d] = sort(f(id,kk));
        accum(id(d(1)),1) = accum(id(d(1)),1) + inf;
        accum(id(d(end)),1) = accum(id(d(end)),1) + inf;
        for jj = 2:length(d)-1
           accum(id(d(jj)),1) =  accum(id(d(jj)),1) + (c(jj+1) - c(jj-1))/(fmax(kk) - fmin(kk)); 
        end
    end
end