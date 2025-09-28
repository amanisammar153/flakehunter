function [mask,labelMap,confMap,d2min, usedParams] = adaptiveInfer(Icorr, P, varargin)
% ADAPTIVEINFER  Runs inferGMM with your defaults, then nudges thresholds
% on outliers only. Returns the final mask + the thresholds actually used.

ip = inputParser;
ip.addParameter('UseRadiusGating', true, @islogical);
ip.addParameter('MaxFGpct', 20, @isscalar);   % if FG > this, tighten
ip.addParameter('MinFGpct', 0.3, @isscalar);  % if 0 < FG < this, relax
ip.addParameter('TightenStep', struct('std',-0.3,'conf',+0.05));
ip.addParameter('RelaxStep',   struct('std',+0.3,'conf',-0.05));
ip.parse(varargin{:});
opt = ip.Results;

% start with project defaults
stdT  = P.stdThresh;
confT = P.confThresh;

% 1st pass
[mask,labelMap,confMap,d2min] = detector.inferGMM(Icorr, P.gmm, ...
    'UsedChannels', P.usedChannels, ...
    'StdThresh', stdT, 'ConfThresh', confT, ...
    'UseRadiusGating', opt.UseRadiusGating, 'ReturnLabels0ForBG', true);

fg = 100*nnz(mask)/numel(mask);

% If way too big: tighten a bit and retry once
if fg > opt.MaxFGpct
    stdT  = max(2.5, stdT + opt.TightenStep.std);
    confT = min(0.95, confT + opt.TightenStep.conf);
    [mask,labelMap,confMap,d2min] = detector.inferGMM(Icorr, P.gmm, ...
        'UsedChannels', P.usedChannels, ...
        'StdThresh', stdT, 'ConfThresh', confT, ...
        'UseRadiusGating', opt.UseRadiusGating, 'ReturnLabels0ForBG', true);
end

% If almost nothing (but not zero): relax a bit and retry once
if fg > 0 && fg < opt.MinFGpct
    stdT  = min(6.0, stdT + opt.RelaxStep.std);
    confT = max(0.40, confT + opt.RelaxStep.conf);
    [mask,labelMap,confMap,d2min] = detector.inferGMM(Icorr, P.gmm, ...
        'UsedChannels', P.usedChannels, ...
        'StdThresh', stdT, 'ConfThresh', confT, ...
        'UseRadiusGating', opt.UseRadiusGating, 'ReturnLabels0ForBG', true);
end

usedParams = struct('StdThresh',stdT,'ConfThresh',confT);
end
