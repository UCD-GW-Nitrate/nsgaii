function opt = nsgaiiOptions(problemSize)

opt.PopulationSize = 50;
opt.ProblemSize = problemSize;
opt.Xmin = zeros(1,problemSize);
opt.Xmax = ones(1,problemSize);
opt.Selection.Method = 'Tournament';
opt.Selection.TourSize = 4;
opt.Crossover.Method = 'Heuristic'; % Uniform Intermediate Heuristic
opt.Crossover.Dinterm = 0.25;
opt.Crossover.Sheuristic = 0.25;