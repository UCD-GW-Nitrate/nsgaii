%% Test Code with test functions
% it's a good practice to clear the work space after running each example.
% Each example can be executed as code block by Ctrl+Enter
clear
clc
%% Schaffer test
nsga_opt = nsgaiiOptions(1);
nsga_opt.Xmin = -1000;
nsga_opt.Xmax = 1000;
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'SCH');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Fonseca test
nsga_opt = nsgaiiOptions(3);
nsga_opt.Xmin = -4*ones(1,3);
nsga_opt.Xmax = 4*ones(1,3);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'FON');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Poloni test
nsga_opt = nsgaiiOptions(2);
nsga_opt.Xmin = -pi*ones(1,2);
nsga_opt.Xmax = pi*ones(1,2);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'POL');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Kursawe test
nsga_opt = nsgaiiOptions(3);
nsga_opt.PopulationSize = 50;
nsga_opt.Xmin = -5*ones(1,3);
nsga_opt.Xmax = 5*ones(1,3);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'KUR');
ParetoSolutions = nsgaii(f, nsga_opt);
%% Zitzler 1 test
nsga_opt = nsgaiiOptions(30);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'ZDT1');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Zitzler 2 test
nsga_opt = nsgaiiOptions(30);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'ZDT2');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Zitzler 3 test
nsga_opt = nsgaiiOptions(30);
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'ZDT3');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Zitzler 4 test
nsga_opt = nsgaiiOptions(10);
nsga_opt.Xmin = [0 -5*ones(1,9)];
nsga_opt.Xmax = [1 5*ones(1,9)];
nsga_opt.MaxGeneration = 1000;
nsga_opt.PlotAllFronts = false;
f = @(x)testObjectiveFunctions(x, 'ZDT4');
ParetoSolutions = nsgaii(f, nsga_opt);

%% Zitzler 6 test
nsga_opt = nsgaiiOptions(10);
nsga_opt.PlotAllFronts = false;
nsga_opt.MaxGeneration = 1000;
f = @(x)testObjectiveFunctions(x, 'ZDT6');
ParetoSolutions = nsgaii(f, nsga_opt);