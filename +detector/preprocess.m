function Icorr = preprocess(Irgb, flatfield, varargin)
% DETECTOR.PREPROCESS  Flat-field / vignette correction (repo-compatible).
%   Icorr = detector.preprocess(Irgb, flatfield, 'FlatfieldMean', [R G B], ...)
%
% Inputs
%   Irgb      : HxWx3 RGB (uint8/uint16/single/double)
%   flatfield : HxWx3 or 1x1x3 RGB flat-field image
%
% Name-Value
%   'FlatfieldMean'      1x3 vector (default: mean over flatfield per channel)
%   'MaxBackgroundValue' scalar clip ceiling (default 241)   % [] to skip
%   'ClipRange'          [lo hi] final clamp (default [0 255])
%   'OutputType'         'uint8' | 'like' | 'double' (default 'uint8')
%   'Epsilon'            denom guard (default 1)             % >0
%   'ROI'                HxW logical mask (optional)         % only correct inside ROI
%
% Output
%   Icorr : corrected image (type per OutputType)

% ---------- parse ----------
ip = inputParser; ip.CaseSensitive = false;
ip.addParameter('FlatfieldMean', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==3));
ip.addParameter('MaxBackgroundValue', 241, @(x) isempty(x) || (isscalar(x)&&isnumeric(x)));
ip.addParameter('ClipRange', [0 255], @(x) isnumeric(x)&&numel(x)==2);
ip.addParameter('OutputType', 'uint8', @(s) ischar(s)||isstring(s));
ip.addParameter('Epsilon', 1, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('ROI', [], @(m) isempty(m) || islogical(m));
ip.parse(varargin{:});
opt = ip.Results;

% ---------- normalize inputs ----------
validateattributes(Irgb, {'numeric'}, {'nonempty','ndims',3});
assert(size(Irgb,3)==3, 'Irgb must be HxWx3.');
H = size(Irgb,1); W = size(Irgb,2);
inClass = class(Irgb);
I = double(Irgb);

% flatfield: allow 1x1x3 or HxWx3
validateattributes(flatfield, {'numeric'}, {'nonempty','ndims',3});
assert(size(flatfield,3)==3, 'flatfield must be HxWx3 or 1x1x3.');
if size(flatfield,1)==1 && size(flatfield,2)==1
    F = repmat(double(flatfield), H, W);
else
    assert(isequal(size(flatfield,1),H) && isequal(size(flatfield,2),W), ...
        'Flatfield size must be 1x1x3 or HxWx3 matching Irgb.');
    F = double(flatfield);
end

% ROI handling
if ~isempty(opt.ROI)
    assert(isequal(size(opt.ROI), [H W]), 'ROI must be HxW logical.');
    ROI = logical(opt.ROI);
else
    ROI = true(H,W);
end

% choose per-channel reference means
if isempty(opt.FlatfieldMean)
    % mean over flatfield under ROI (robust to zeros / mask)
    mu = zeros(1,3);
    for c=1:3
        v = F(:,:,c);
        vc = v(ROI);
        vc = vc(isfinite(vc));
        if isempty(vc), vc = v(:); end
        mu(c) = mean(vc(:), 'omitnan');
    end
else
    mu = double(opt.FlatfieldMean(:)).';
end

epsd = max(1e-12, double(opt.Epsilon));

% ---------- correction ----------
Icorr = zeros(H,W,3);
for c = 1:3
    denom = max(F(:,:,c), epsd);    % guard division
    val = I(:,:,c) ./ denom * mu(c);
    % optional background ceiling (clip) INSIDE ROI only
    if ~isempty(opt.MaxBackgroundValue)
        val(ROI) = min(val(ROI), double(opt.MaxBackgroundValue));
    end
    Icorr(:,:,c) = val;
end

% final clamp and cast
lo = opt.ClipRange(1); hi = opt.ClipRange(2);
Icorr = max(min(Icorr, hi), lo);

switch lower(string(opt.OutputType))
    case "like",  Icorr = cast(Icorr, inClass);
    case "uint8", Icorr = uint8(Icorr);
    case "double" % leave as double
    otherwise, error('Unsupported OutputType: %s', opt.OutputType);
end
end
