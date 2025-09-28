function PR = eval_computePR(scoresPos, scoresNeg, varargin)
% DETECTOR.EVAL_COMPUTEPR  Precision/Recall/F1 vs threshold from scores.
%
% PR = detector.eval_computePR(scoresPos, scoresNeg, 'Thresh', linspace(0,1,101))
%
% Inputs
%   scoresPos : vector of scores for positives   (e.g., GT pixels=1 or GT flakes)
%   scoresNeg : vector of scores for negatives   (e.g., GT pixels=0)
% Name–Value
%   'Thresh'  : thresholds to evaluate (default 0:0.01:1)
%
% Output struct PR with fields:
%   .thr, .tp, .fp, .fn, .tn, .prec, .rec, .f1, .ap
%
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Thresh', linspace(0,1,101), @(x)isnumeric(x)&&isvector(x));
ip.parse(varargin{:});
thr = sort(ip.Results.Thresh(:).');

pos = double(scoresPos(:));
neg = double(scoresNeg(:));
pos(isnan(pos)) = []; neg(isnan(neg)) = [];

np = numel(pos); nn = numel(neg);
tp = zeros(size(thr)); fp = zeros(size(thr));
for i = 1:numel(thr)
    t = thr(i);
    tp(i) = sum(pos >= t);
    fp(i) = sum(neg >= t);
end
fn = np - tp;
tn = nn - fp;

prec = tp ./ max(1, tp+fp);
rec  = tp ./ max(1, tp+fn);
f1   = 2*(prec.*rec) ./ max(1e-12, prec+rec);

% Average Precision via step-wise integration over recall
[recS, idx] = sort(rec, 'ascend'); precS = prec(idx);
ap = trapz(recS, precS);

PR = struct('thr',thr,'tp',tp,'fp',fp,'fn',fn,'tn',tn, ...
            'prec',prec,'rec',rec,'f1',f1,'ap',ap);
end
