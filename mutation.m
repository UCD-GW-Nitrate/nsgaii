function xnew = mutation(x, options)
%disp('Mutating...');
xnew = x;
switch options.Mutation.Method
    case 'Uniform'
        for ii = 1:size(x,2)
            test = rand(options.PopulationSize,1) < options.Mutation.Probability; 
            xnew(test,ii) = options.Xmin(ii) + (options.Xmax(ii) - options.Xmin(ii))*rand(sum(test),1);
            test = xnew(:,ii) < options.Xmin(ii);
            xnew(test,ii) = options.Xmin(ii);
            test = xnew(:,ii) > options.Xmax(ii);
            xnew(test,ii) = options.Xmax(ii);
        end
    case 'Gaussian'
        warning('Not debuged yet')
        options.Mutation.ScaleGaussian = options.Mutation.ScaleGaussian * ...
            (1-options.Mutation.ShrinkGaussian*options.Generation/options.MaxGeneration);
        for ii = 1:size(x,2)
            test = rand(options.PopulationSize,1) < options.Mutation.Probability;
             xnew(test,ii) =  xnew(test,ii) + ...
                 normrnd(0,options.Mutation.ScaleGaussian, sum(test), 1);
        end
        
    case 'Breeder'
        warning('Not debuged yet')
        r = (options.Generation - 1)*(0.05 - 0.5)/(options.MaxGeneration - 1) + 0.5;
        k = (options.Generation - 1)*(8 - 2)/(options.MaxGeneration - 1) + 2;
        for ii = 1:size(x,2)
            test = rand(options.PopulationSize,1) < options.Mutation.Probability;
            nMutations = sum(test);
            if nMutations > 0
                a = 2^(-rand(nMutations,1)*k);
                s = -1+2*rand(nMutations,1);
                ri = r*(options.Xmax - options.Xmin);
                xnew(test,ii) =  xnew(test,ii) + s.*ri.*a;
            end
            test = xnew(:,ii) < options.Xmin(ii);
            xnew(test,ii) = options.Xmin(ii);
            test = xnew(:,ii) > options.Xmax(ii);
            xnew(test,ii) = options.Xmax(ii);
            
        end
end


