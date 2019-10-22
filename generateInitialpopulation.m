function Population = generateInitialpopulation(options)

Xmin = repmat(options.Xmin, options.PopulationSize,1);
Xmax = repmat(options.Xmax, options.PopulationSize,1);

Population = Xmin + (Xmax - Xmin).*rand(options.PopulationSize, options.ProblemSize);