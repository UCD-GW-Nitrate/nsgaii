function idSelected = selectSurvivors(rank, accum, options)
%disp('Selecting survivors...');

idSelected = [];
avaliableRanks = unique(rank);
for ii = 1:length(avaliableRanks)
    id = find(rank == avaliableRanks(ii));
    if length(idSelected) + length(id) <= options.PopulationSize
        idSelected = [idSelected; id];
    else
        nSelected = options.PopulationSize - length(idSelected);
        [~, d] = sort(accum(id),'descend');
        idSelected = [idSelected; id(d(1:nSelected))];       
    end
end