function REP = eval_regionMetrics(detMasks, gtMasks, varargin)
% DETECTOR.EVAL_REGIONMETRICS  Aggregate region-level metrics across images.
%
% REP = detector.eval_regionMetrics(detMasks, gtMasks, ...
%           'IoUThresh',0.5, 'PerImage',true)
%
% Inputs
%   detMasks, gtMasks : HxW logical or cell arrays (same count)
%
% Name–Value
%   'IoUThresh' : scalar (default 0.5)
%   'PerImage'  : logical — include per-image table (default true)
%
% Output
%   REP.micro : struct with TP,FP,FN, precision, recall, f1 (micro-averaged)
%   REP.macro : struct with mean precision/recall/f1 across images
%   REP.perImage : table with metrics per image (if requested)
%
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('IoUThresh',0.5,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('PerImage',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.parse(varargin{:});
opt = ip.Results;

D = detMasks; G = gtMasks;
if ~iscell(D), D={D}; end
if ~iscell(G), G={G}; end
assert(numel(D)==numel(G),'Counts must match.');

rows = [];
TP=0; FP=0; FN=0;
for i=1:numel(D)
    Mi = detector.eval_matchRegions(D{i}, G{i}, 'IoUThresh', opt.IoUThresh);
    TP = TP + Mi.TP; FP = FP + Mi.FP; FN = FN + Mi.FN;
    rows = [rows; i, Mi.TP, Mi.FP, Mi.FN, Mi.precision, Mi.recall, Mi.f1]; %#ok<AGROW>
end

prec_micro = TP / max(1,TP+FP);
rec_micro  = TP / max(1,TP+FN);
f1_micro   = 2*prec_micro*rec_micro / max(1e-12, prec_micro+rec_micro);

T = array2table(rows, 'VariableNames', {'image','TP','FP','FN','precision','recall','f1'});

REP.micro = struct('TP',TP,'FP',FP,'FN',FN,'precision',prec_micro,'recall',rec_micro,'f1',f1_micro);
REP.macro = struct('precision', mean(T.precision), 'recall', mean(T.recall), 'f1', mean(T.f1));
REP.perImage = tern(opt.PerImage, T, table());
end

function x = tern(c,a,b), if c, x=a; else, x=b; end
end
