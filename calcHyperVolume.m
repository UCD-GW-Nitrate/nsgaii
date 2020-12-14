function hv = calcHyperVolume(obj, far_p, uto_p)
% hv = calcHyperVolume(obj, far_p)
% Calculates the hypervolume of the pareto front described by the variable
% obj. The goal for the objective functions is assumed to be the
% minimization. If not negate the ones that are maximized
%
% obj is a matrix with dimensions Nsol x Nobj, where Nsol is the number of
% solutions and Nobj is the number of objectives. Currently only 2
% objective pareto front is implemented
% far_p is the point that is dominated by every point of the pareto front.
% uto_p is the point that dominates every point of the pareto front.


if size(obj,2) ~= 2
    error('This function is not defined for other than 2D pareto fronts');
end

% Normalize the objective function between 0 and 1 based on the far
obj1 = (obj(:,1) - uto_p(1))/(far_p(1) - uto_p(1));
obj2 = (obj(:,2) - uto_p(2))/(far_p(2) - uto_p(2));

% Sort the pareto front based on the second objective function, eg the y
[obj2, ind] = sort(obj2, 'descend');
obj1 = obj1(ind,1);
hv = 0;
for ii = 1:length(obj2)
    if ii == 1
        hv = (1 - obj1(ii))*(1 - obj2(ii));
    else
        hv = hv+ (1 - obj1(ii))*(obj2(ii-1) - obj2(ii));
    end
end
