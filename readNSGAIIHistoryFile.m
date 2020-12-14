function ParetoHistory = readNSGAIIHistoryFile(varargin)
if length(varargin) == 1
    filename = varargin{1};
    nalloc = 10000;
end
if length(varargin) == 2 
    filename = varargin{1};
    nalloc = varargin{2};
end

fid = fopen(filename,'r');

nobj = 0;
ParetoHistory.ObjLim = zeros(2);
ParetoHistory.Nsol = zeros(nalloc,1);
ParetoHistory.Gen(nalloc,1).Obj = [];
cnt_gen = 1;
while ~feof(fid)
    tline = fgetl(fid);
    if isempty(tline)
        break;
    end
    gen_info = textscan(tline, '%f',3);
    nobj = gen_info{1}(2);
    nsol = gen_info{1}(3);
    ParetoHistory.Nsol(cnt_gen,1) = nsol;
    sol = textscan(fid, '%f',nobj*nsol);
    sol = reshape(sol{1},nobj, nsol)';
    ParetoHistory.Gen(cnt_gen,1).Obj = sol;
    if cnt_gen == 1
        ParetoHistory.ObjLim(1,:) = min(sol,[],1);
        ParetoHistory.ObjLim(2,:) = max(sol,[],1);
    else
        ParetoHistory.ObjLim(1,:) = min(ParetoHistory.ObjLim(1,:), min(sol,[],1));
        ParetoHistory.ObjLim(2,:) = max(ParetoHistory.ObjLim(2,:), max(sol,[],1));
    end
    cnt_gen = cnt_gen + 1;
    tline = fgetl(fid);
end

fclose(fid);
ParetoHistory.Gen(cnt_gen:end,:) = [];
ParetoHistory.Nsol(cnt_gen:end,:) = [];