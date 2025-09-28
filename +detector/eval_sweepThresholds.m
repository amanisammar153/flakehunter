function REP = eval_sweepThresholds(confMaps, gtMasks, varargin)
% DETECTOR.EVAL_SWEEPTHRESHOLDS  Sweep ConfThresh across images and report PR.
%
% REP = detector.eval_sweepThresholds(confMaps, gtMasks, ...
%          'Thresh', 0:0.01:1, 'ROIs',[], 'SaveCSV','', 'Plot',false)
%
% Inputs
%   confMaps : HxW numeric, or cell array of HxW maps in [0,1]
%   gtMasks  : HxW logical, or cell array of HxW logical masks (same size as confMaps)
%
% Name–Value
%   'Thresh'  : thresholds to evaluate (default 0:0.01:1)
%   'ROIs'    : [] or single HxW logical or cell array (same length as confMaps)
%   'SaveCSV' : '' or path to write a CSV with PR curve (default '')
%   'Plot'    : logical to show a PR curve (default false)
%
% Output struct REP:
%   .PR            : struct from eval_computePR (thr, prec, rec, f1, ap, ...)
%   .bestIdx       : index into thr of best F1
%   .bestThreshold : threshold value
%   .bestF1        : best F1
%   .bestPrec      : precision at best
%   .bestRec       : recall at best
%
% Example
%   REP = detector.eval_sweepThresholds({conf1,conf2},{gt1,gt2}, 'Thresh',0:0.02:1,'SaveCSV','pr.csv');

ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Thresh', linspace(0,1,101), @(x)isnumeric(x)&&isvector(x));
ip.addParameter('ROIs', [], @(r) isempty(r) || islogical(r) || iscell(r));
ip.addParameter('SaveCSV','', @(s)ischar(s)||isstring(s));
ip.addParameter('Plot', false, @(b)islogical(b)||ismember(b,[0 1]));
ip.parse(varargin{:});
opt = ip.Results;

% normalize inputs to cell arrays
C = confMaps; G = gtMasks; R = opt.ROIs;
if ~iscell(C), C = {C}; end
if ~iscell(G), G = {G}; end
assert(numel(C)==numel(G), 'confMaps and gtMasks must have same count.');
if isempty(R)
    R = cell(size(C));  % all empty
elseif ~iscell(R)
    R = repmat({R}, size(C));
else
    assert(numel(R)==numel(C), 'ROIs length must match confMaps.');
end

% aggregate scores
sp_all = []; sn_all = [];
for i=1:numel(C)
    [sp, sn] = detector.eval_masks2scores(C{i}, G{i}, 'ROI', R{i});
    sp_all = [sp_all; sp]; %#ok<AGROW>
    sn_all = [sn_all; sn]; %#ok<AGROW>
end

% compute PR across thresholds
PR = detector.eval_computePR(sp_all, sn_all, 'Thresh', opt.Thresh);

% pick best F1 (ties → lowest threshold)
[bestF1, bi] = max(PR.f1);
bestThr = PR.thr(bi);
bestPrec = PR.prec(bi);
bestRec  = PR.rec(bi);

% optional CSV
if ~isempty(opt.SaveCSV)
    T = table(PR.thr(:), PR.prec(:), PR.rec(:), PR.f1(:), ...
        'VariableNames', {'threshold','precision','recall','f1'});
    try
        writetable(T, char(opt.SaveCSV));
    catch ME
        warning('eval_sweepThresholds:CSV','Failed to save CSV (%s): %s', opt.SaveCSV, ME.message);
    end
end

% optional plot
if opt.Plot
    figure('Name','PR curve','NumberTitle','off');
    plot(PR.rec, PR.prec, '-', 'LineWidth',1.5); grid on;
    xlabel('Recall'); ylabel('Precision'); title(sprintf('PR (AP=%.3f)', PR.ap));
    hold on; plot(bestRec, bestPrec, 'o', 'MarkerSize',6, 'LineWidth',1.5);
end

% pack
REP = struct('PR',PR, 'bestIdx',bi, 'bestThreshold',bestThr, ...
             'bestF1',bestF1, 'bestPrec',bestPrec, 'bestRec',bestRec);
end
