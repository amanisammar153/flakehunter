function GS = eval_gridSearch(confMaps, gtMasks, varargin)
% DETECTOR.EVAL_GRIDSEARCH  Tiny grid-search over postprocess knobs.
%
% GS = detector.eval_gridSearch(confMaps, gtMasks, ...
%         'ThreshList', 0:0.05:1, 'MinArea_um2_List', [0 5 10], ...
%         'UmPerPixel', 1, 'OpenRadiusPx',0, 'CloseRadiusPx',0, ...
%         'MinSolidity',0, 'KeepTopK',0)
%
% Inputs
%   confMaps, gtMasks : HxW arrays or cell arrays of equal length
%
% Name–Value (all optional)
%   'ThreshList'       : vector of thresholds (default 0:0.05:1)
%   'MinArea_um2_List' : vector of min areas in µm² (default [0])
%   'UmPerPixel'       : scalar µm/px (default 1)   % used to convert µm² → px
%   'OpenRadiusPx'     : int (default 0)
%   'CloseRadiusPx'    : int (default 0)
%   'MinSolidity'      : [0..1] (default 0)
%   'KeepTopK'         : integer (default 0)
%
% Output struct GS
%   .table   : leaderboard table sorted by F1 desc
%   .best    : struct with best params and metrics (F1, precision, recall)
%   .details : per-combo pixel counts (TP/FP/FN/TN)
%
% Notes
% - Metrics are pixel-level: we form a mask by postprocessing (with the
%   given knobs) and compare to gt mask as positives.
% - If you want region-level metrics, we can add an overlap-based variant next.

% ----- parse -----
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('ThreshList', 0:0.05:1, @(x)isnumeric(x)&&isvector(x));
ip.addParameter('MinArea_um2_List', 0, @(x)isnumeric(x)&&isvector(x));
ip.addParameter('UmPerPixel', 1, @(x)isscalar(x)&&x>0);
ip.addParameter('OpenRadiusPx', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('CloseRadiusPx', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('MinSolidity', 0, @(x)isscalar(x)&&x>=0);
ip.addParameter('KeepTopK', 0, @(x)isscalar(x)&&x>=0);
ip.parse(varargin{:});
opt = ip.Results;

% normalize to cell
C = confMaps; G = gtMasks;
if ~iscell(C), C={C}; end
if ~iscell(G), G={G}; end
assert(numel(C)==numel(G), 'confMaps and gtMasks must have same count.');

% precompute pixels total within evaluation area
Nimgs = numel(C);
results = [];  % will collect rows

for a = 1:numel(opt.MinArea_um2_List)
    minA = opt.MinArea_um2_List(a);
    for t = 1:numel(opt.ThreshList)
        thr = opt.ThreshList(t);

        TP=0; FP=0; FN=0; TN=0;
        for i=1:Nimgs
            conf = single(C{i});
            gt   = logical(G{i});
            H = size(conf,1); W = size(conf,2);
            labels = zeros(H,W,'uint16'); % not used, but postprocess expects it

            % postprocess with current knobs
            pp = detector.postprocess(conf, labels, struct( ...
                'ConfThresh', thr, ...
                'OpenRadiusPx', opt.OpenRadiusPx, ...
                'CloseRadiusPx', opt.CloseRadiusPx, ...
                'MinArea_um2', minA, ...
                'AutoAreaPct', [], ...
                'umPerPixel', opt.UmPerPixel, ...
                'MinSolidity', opt.MinSolidity, ...
                'KeepTopK', opt.KeepTopK ));

            pred = logical(pp.mask);

            % pixel confusion
            TP = TP + sum(pred(:) &  gt(:));
            FP = FP + sum(pred(:) & ~gt(:));
            FN = FN + sum(~pred(:) &  gt(:));
            TN = TN + sum(~pred(:) & ~gt(:));
        end

        prec = TP / max(1, TP+FP);
        rec  = TP / max(1, TP+FN);
        f1   = 2*(prec*rec)/max(1e-12, (prec+rec));

        results = [results; thr, minA, prec, rec, f1, TP, FP, FN, TN]; %#ok<AGROW>
    end
end

T = array2table(results, 'VariableNames', ...
    {'threshold','min_area_um2','precision','recall','f1','TP','FP','FN','TN'});
% sort by F1 desc, tie-break: higher precision, then lower threshold
T = sortrows(T, {'f1','precision','threshold'}, {'descend','descend','ascend'});

best = table2struct(T(1,:));

GS = struct('table', T, 'best', best, ...
    'details', struct('TP',best.TP,'FP',best.FP,'FN',best.FN,'TN',best.TN));
end
