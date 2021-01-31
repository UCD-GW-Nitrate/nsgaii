function ParetoSolutions = readNSGAIIOutputFile(filename)
ParetoSolutions.x = [];
ParetoSolutions.f = [];
fid = fopen(filename,'r');
C = textscan(fid, '%f',3);
Nsol = C{1}(1);
Nvar = C{1}(2);
Nobj = C{1}(3);
ParetoSolutions.x = nan(Nsol, Nvar);
ParetoSolutions.f = nan(Nsol, Nobj);
% read Desicion variables 
for ii = 1:Nsol
    tmp = textscan(fid, '%f',Nvar);
    ParetoSolutions.x(ii,:) = tmp{1,1}';
end

for ii = 1:Nsol
    tmp = textscan(fid, '%f',Nobj);
    ParetoSolutions.f(ii,:) = tmp{1,1}';
end
fclose(fid);