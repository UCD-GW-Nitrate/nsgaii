function xnew = crossoverParents(x, id, rank, accum, options)
%disp('Crossover...');

switch options.Crossover.Method
    case 'Uniform'
        for ii = 1:2:options.PopulationSize
            if ii+1 > options.PopulationSize
                break;
            end
            test = rand(1, options.ProblemSize) > 0.5;
            xnew(ii, test) = x(id(ii), test);
            xnew(ii, ~test) = x(id(ii+1), ~test);
            xnew(ii+1, test) = x(id(ii+1), test);
            xnew(ii+1, ~test) = x(id(ii), ~test);
        end
    case 'Intermediate'
        for ii = 1:2:options.PopulationSize
            if ii+1 > options.PopulationSize
                break;
            end
            r = -opt.Crossover.Dinterm + (1+2*opt.Crossover.Dinterm)*rand(1, options.ProblemSize);
            xnew(ii,:) = r*x(id(ii),:) + (1-r)*x(id(ii+1),:);
            xnew(ii+1,:) = (1-r)*x(id(ii),:) + r*x(id(ii+1),:);
        end
    case 'Heuristic'
        ii = 1;
        for jj = 1:2:length(id)
            if rank(id(jj)) < rank(id(jj+1))
                p1 = x(id(jj),:);
                p2 = x(id(jj+1),:);
            elseif rank(id(jj)) > rank(id(jj+1))
                p2 = x(id(jj),:);
                p1 = x(id(jj+1),:);
            elseif rank(id(jj)) == rank(id(jj+1))
                if accum(id(jj)) >= accum(id(jj+1))
                    p1 = x(id(jj),:);
                    p2 = x(id(jj+1),:);
                else
                    p2 = x(id(jj),:);
                    p1 = x(id(jj+1),:);
                end
            end
            
            r = abs(normrnd(0,options.Crossover.Sheuristic,1,options.ProblemSize));
            xnew(ii,:) = p1 + r.*(p1-p2);
            for k = 1:options.ProblemSize
                if xnew(ii,k) > options.Xmin(k) && xnew(ii,k) < options.Xmax(k)
                    continue;
                end
                s = options.Crossover.Sheuristic;
                while true
                    r = abs(normrnd(0,s));
                    if s < 0
                        if rand > 0.5
                            r = 0;
                        else
                            if newvalue < options.Xmin(k)
                                 xnew(ii,k) = options.Xmin(k);
                            elseif newvalue > options.Xmax(k)
                                xnew(ii,k) = options.Xmax(k);
                            end
                            break;
                        end
                            
                    end
                    newvalue = p1(k) + r*(p1(k)-p2(k));
                    if newvalue > options.Xmin(k) && newvalue < options.Xmax(k)
                        xnew(ii,k) = newvalue;
                        break;
                    end
                    s = s-0.01;
                end
            end
                
            
            %test = xnew(ii,:) > options.Xmax;
            %xnew(ii,test) = options.Xmax(test);
            %test = xnew(ii,:) < options.Xmin;
            %xnew(ii,test) = options.Xmin(test);
            ii = ii + 1;
        end
        
end
