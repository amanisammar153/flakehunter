function OUT = eval_datasetTuner(confMaps, gtMasks, outDir, varargin)
% DETECTOR.EVAL_DATASETTUNER  Sweep thresholds, rank combos, and save report.
%
% OUT = detector.eval_datasetTuner(confMaps, gtMasks, outDir, ...
%       'Thresh',0:0.02:1, 'MinArea_um2_List',[0 50 200], 'UmPerPixel',1, ...
%       'OpenRadiusPx',0, 'CloseRadiusPx',0, 'MinSolidity',0, 'KeepTopK',0, ...
%       'Pick','bestF1', ...             % 'bestF1'|'targetPrec'|'targetRec'
%       'Target',0.9, ...                % for target modes
%       'SavePrefix','tune')
%
% INPUTS (arrays or cell arrays)
%   confMaps : HxW numeric in [0,1] or cell of them
%   gtMasks  : HxW logical          or cell of them
%   outDir   : folder to write CSV/PNG/JSON
%
% OUTPUT
%   OUT.tableCSV   : path to leaderboard CSV
%   OUT.prCSV      : path to PR CSV (all images pooled)
%   OUT.prPNG      : path to PR PNG
%   OUT.pick       : struct with chosen settings (threshold, min_area_um2, etc.)
%
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Thresh', 0:0.02:1, @(x)isnumeric(x)&&isvector(x));
ip.addParameter('MinArea_um2_List', 0, @(x)isnumeric(x)&&isvector(x));
ip.addParameter('UmPerPixel', 1, @(x)isscalar(x)&&x>0);
ip.addParameter('OpenRadiusPx', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('CloseRadiusPx', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('MinSolidity', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('KeepTopK', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('Pick','bestF1',@(s)any(strcmpi(s,{'bestF1','targetPrec','targetRec'})));
ip.addParameter('Target', 0.9, @(x)isscalar(x)&&x>=0&&x<=1);
ip.addParameter('SavePrefix','tune',@(s)ischar(s)||isstring(s));
ip.parse(varargin{:});
opt = ip.Results;

if ~isfolder(outDir), mkdir(outDir); end
pref = char(opt.SavePrefix);

% 1) build pooled PR (pixel-level) for plotting/report
C = confMaps; G = gtMasks;
if ~iscell(C), C={C}; end
if ~iscell(G), G={G}; end
sp_all = []; sn_all = [];
for i=1:numel(C)
    [sp,sn] = detector.eval_masks2scores(C{i}, logical(G{i}));
    sp_all = [sp_all; sp]; %#ok<AGROW>
    sn_all = [sn_all; sn]; %#ok<AGROW>
end
PR = detector.eval_computePR(sp_all, sn_all, 'Thresh', opt.Thresh);

% Save PR CSV/PNG
prCSV = fullfile(outDir, sprintf('%s_pr.csv', pref));
writetable(table(PR.thr(:),PR.prec(:),PR.rec(:),PR.f1(:), ...
    'VariableNames',{'threshold','precision','recall','f1'}), prCSV);

prPNG = fullfile(outDir, sprintf('%s_pr.png', pref));
figure('Visible','off'); plot(PR.rec, PR.prec, '-', 'LineWidth',1.5); grid on;
xlabel('Recall'); ylabel('Precision'); title(sprintf('Pooled PR (AP=%.3f)', PR.ap));
saveas(gcf, prPNG); close(gcf);

% 2) grid search over postprocess knobs
GS = detector.eval_gridSearch(C, G, ...
    'ThreshList', opt.Thresh, ...
    'MinArea_um2_List', opt.MinArea_um2_List, ...
    'UmPerPixel', opt.UmPerPixel, ...
    'OpenRadiusPx', opt.OpenRadiusPx, ...
    'CloseRadiusPx', opt.CloseRadiusPx, ...
    'MinSolidity', opt.MinSolidity, ...
    'KeepTopK', opt.KeepTopK);

tableCSV = fullfile(outDir, sprintf('%s_leaderboard.csv', pref));
writetable(GS.table, tableCSV);

% 3) choose settings
switch lower(opt.Pick)
    case 'bestf1'
        pickRow = GS.table(1,:);
    case 'targetprec'
        % highest recall with precision >= Target; tie-break lower threshold
        T = GS.table(GS.table.precision >= opt.Target, :);
        if isempty(T), T = GS.table; end
        T = sortrows(T, {'recall','threshold'},{'descend','ascend'});
        pickRow = T(1,:);
    case 'targetrec'
        % highest precision with recall >= Target; tie-break lower threshold
        T = GS.table(GS.table.recall >= opt.Target, :);
        if isempty(T), T = GS.table; end
        T = sortrows(T, {'precision','threshold'},{'descend','ascend'});
        pickRow = T(1,:);
end

pick = table2struct(pickRow);
pick.AP_pooled = PR.ap;

% 4) save pick JSON
pickJSON = fullfile(outDir, sprintf('%s_pick.json', pref));
fid=fopen(pickJSON,'w'); fwrite(fid, jsonencode(pick),'char'); fclose(fid);

OUT = struct('tableCSV',tableCSV,'prCSV',prCSV,'prPNG',prPNG,'pick',pick,'pickJSON',pickJSON);
end
