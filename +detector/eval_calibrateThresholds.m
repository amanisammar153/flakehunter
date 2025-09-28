function pick = eval_calibrateThresholds(scoresPos, scoresNeg, varargin)
% DETECTOR.EVAL_CALIBRATETHRESHOLDS  Pick threshold for target prec/rec or best-F1.
%
% pick = detector.eval_calibrateThresholds(sp, sn, ...
%          'Mode','bestF1'|'targetPrec'|'targetRec', 'Target',0.9, 'Thresh',0:0.005:1)
%
% Output 'pick' struct: threshold, precision, recall, f1, ap, idx
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Mode','bestF1',@(s)any(strcmpi(s,{'bestF1','targetPrec','targetRec'})));
ip.addParameter('Target',0.9,@(x)isscalar(x)&&x>=0&&x<=1);
ip.addParameter('Thresh', 0:0.005:1, @(x)isnumeric(x)&&isvector(x));
ip.parse(varargin{:});
opt = ip.Results;

PR = detector.eval_computePR(scoresPos, scoresNeg, 'Thresh', opt.Thresh);

switch lower(opt.Mode)
    case 'bestf1'
        [~,idx] = max(PR.f1);
    case 'targetprec'
        idx = find(PR.prec >= opt.Target, 1, 'last'); % prefer higher threshold
        if isempty(idx), [~,idx] = max(PR.prec); end
    case 'targetrec'
        idx = find(PR.rec >= opt.Target, 1, 'first'); % prefer lower threshold
        if isempty(idx), [~,idx] = max(PR.rec); end
end

pick = struct('threshold', PR.thr(idx), 'precision', PR.prec(idx), ...
              'recall', PR.rec(idx), 'f1', PR.f1(idx), 'ap', PR.ap, 'idx', idx);
end
