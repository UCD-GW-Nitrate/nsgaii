function [f, x] = removeDuplicateSolutions(f, x, options)


ii = 1;
while true
    nf = size(f,1);
    if ii == 1
        testSet = 2:nf;
    elseif ii == nf
        testSet = 1:nf-1;
    else
       testSet = [1:ii-1 ii+1:nf]';
    end
    
    dst = sqrt(sum(bsxfun(@minus, f(testSet,:),f(ii,:)).^2, 2));
    id = find(dst < options.UniqueSolutionThreshold);
    if ~isempty(id)
        f(testSet(id),:) = [];
        x(testSet(id),:) = [];
    end
    ii = ii + 1;
    if ii > size(f,1)
        break;
    end
end