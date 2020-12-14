function ParetoSolutions = readNSGAIIOutputFile(filename)
ParetoSolutions.x = [];
ParetoSolutions.f = [];
fid = fopen(filename,'r');
C = textscan(fid, '%d %d %d\n');
frmVar = [];
for ii = 1:double(C{2})
    frmVar = [frmVar '%f '];
end
frmVar = [frmVar '\n'];

frmObj = [];
for ii = 1:double(C{3})
    frmObj = [frmObj '%f '];
end
frmObj = [frmObj '\n'];

for ii = 1:double(C{1})
    V = cell2mat(textscan(fid, frmVar,  1));
    ParetoSolutions.x = [ParetoSolutions.x; V]; 
end

for ii = 1:double(C{1})
    V = cell2mat(textscan(fid, frmObj,  1));
    ParetoSolutions.f = [ParetoSolutions.f; V]; 
end
fclose(fid);