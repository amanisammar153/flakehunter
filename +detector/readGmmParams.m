function out = readGmmParams(src, varargin)
% DETECTOR.READGMMPARAMS  Load & normalize GMM parameters to MATLAB RGB unit scale.
%
% out = detector.readGmmParams(src, 'UsedChannels','RGB', 'AssumeOrder','auto'|'RGB'|'BGR')
%
% INPUT
%   src : filename (.json/.mat) OR struct already in memory
%
% OPTIONS (name–value)
%   'UsedChannels'  char subset of 'RGB' (default 'RGB')
%   'AssumeOrder'   'auto'|'RGB'|'BGR'  (default 'auto')
%                   When 'auto', tries to infer from common Python fields.
%   'PreferInv'     logical, prefer using invSigma if available (default true)
%
% OUTPUT (all single precision, unit scale)
%   out.mu      C x K   (C=#used channels, K=components)          in [0,1]
%   out.Sigma   C x C x K
%   out.radius  C x K   (empty if not provided)                    in [0,1]
%   out.usedStr char like 'RGB','RG','B'
%   out.scaleIn {'unit','uint8'}  original scale that was detected
%   out.channelOrderIn {'RGB','BGR'}  original order detected
%
% Notes
% - If src has *_rgb fields as Kx3, they are transposed internally to 3xK.
% - Covariance scaling: if original scale is uint8 (0..255), Sigma is divided by 255^2.

% ---- parse options ----
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('UsedChannels','RGB',@(s)ischar(s)||isstring(s));
ip.addParameter('AssumeOrder','auto',@(s) any(strcmpi(s,{'auto','RGB','BGR'})));
ip.addParameter('PreferInv',true,@(b)islogical(b)||isscalar(b));
ip.parse(varargin{:});
opt = ip.Results;

% ---- load input struct ----
if ischar(src) || isstring(src)
    fn = char(src);
    assert(exist(fn,'file')==2, 'File not found: %s', fn);
    [~,~,ext] = fileparts(fn);
    switch lower(ext)
        case '.json'
            S = jsondecode(fileread(fn));
        case '.mat'
            tmp = load(fn);
            % take the first struct-like entry
            k = fieldnames(tmp);
            if isempty(k), error('MAT file had no variables.'); end
            S = tmp.(k{1});
        otherwise
            error('Unsupported file extension: %s', ext);
    end
elseif isstruct(src)
    S = src;
else
    error('Unsupported src type. Provide .json/.mat or a struct.');
end

% ---- discover possible field names ----
% Means
if isfield(S,'mu')
    mu = S.mu;                    % expect C x K
elseif isfield(S,'mu_rgb')
    mu = S.mu_rgb.';              % K x 3 -> 3 x K
else
    error('Missing mu / mu_rgb in params.');
end
% Covariances or inverses
invSigma = [];
Sigma    = [];
if isfield(S,'invSigma')
    invSigma = S.invSigma;        % C x C x K
elseif isfield(S,'invSigma_rgb')
    invSigma = S.invSigma_rgb;    % assume already 3x3xK
end
if isempty(invSigma)
    if isfield(S,'Sigma')
        Sigma = S.Sigma;          % C x C x K
    elseif isfield(S,'Sigma_rgb')
        Sigma = S.Sigma_rgb;      % 3 x 3 x K
    else
        error('Missing Sigma / Sigma_rgb / invSigma in params.');
    end
end
% Radii
radius = [];
if isfield(S,'radius')
    radius = S.radius;            % C x K
elseif isfield(S,'radius_rgb')
    radius = S.radius_rgb.';      % K x 3 -> 3 x K
end

% Coerce to double for math
mu = double(mu);
if ~isempty(Sigma), Sigma = double(Sigma); end
if ~isempty(invSigma), invSigma = double(invSigma); end
if ~isempty(radius), radius = double(radius); end

% ---- decide original channel order (RGB vs BGR) ----
orderIn = upper(string(opt.AssumeOrder));
if orderIn=="AUTO"
    % Heuristic: Python exports often use "*_bgr" field names or note "bgr" in comments.
    if any(strcmpi(fieldnames(S),'mu_bgr')) || any(strcmpi(fieldnames(S),'Sigma_bgr')) || any(strcmpi(fieldnames(S),'radius_bgr'))
        orderIn = "BGR";
    else
        % Fall back to RGB (most of your MATLAB inputs)
        orderIn = "RGB";
    end
end

% Reorder to MATLAB RGB if needed
if orderIn=="BGR"
    perm = [3 2 1];   % BGR -> RGB
else
    perm = [1 2 3];
end

% mu: C x K (maybe not 3 if user passed grayscale)
if size(mu,1)==3
    mu = mu(perm, :);
end
% Sigma / invSigma: C x C x K
if ~isempty(Sigma) && size(Sigma,1)==3
    Sigma = Sigma(perm,perm,:);
end
if ~isempty(invSigma) && size(invSigma,1)==3
    invSigma = invSigma(perm,perm,:);
end
% radius: C x K
if ~isempty(radius) && size(radius,1)==3
    radius = radius(perm,:);
end

% ---- detect scale (0..255 vs 0..1) ----
% If any mean > 1.2 or any diag(Sigma) > 1.5, assume uint8 scale.
maxMu   = max(abs(mu(:)));
diagSig = [];
if ~isempty(Sigma)
    for k=1:size(Sigma,3), diagSig(end+1,:) = diag(Sigma(:,:,k)).'; end %#ok<AGROW>
    maxDiag = max(abs(diagSig(:)));
else
    maxDiag = 0;
end
isUint8Scale = (maxMu > 1.2) || (maxDiag > 1.5);  % loose but reliable

scaleIn = tern(isUint8Scale, 'uint8', 'unit');

% Rescale to unit if needed
if isUint8Scale
    mu = mu / 255;
    if ~isempty(Sigma),  Sigma = Sigma / (255^2); end
    if ~isempty(invSigma), invSigma = invSigma * (255^2); end % inv scales opposite
    if ~isempty(radius), radius = radius / 255; end
end

% ---- select UsedChannels subset ----
used = upper(char(opt.UsedChannels));
used = used(ismember(used,'RGB'));
if isempty(used), used = 'RGB'; end
map = struct('R',1,'G',2,'B',3);
idx = arrayfun(@(ch) map.(char(ch)), num2cell(used));

Cfull = size(mu,1);
assert(Cfull>=numel(idx), 'Params have %d channels; UsedChannels asks for %d.', Cfull, numel(idx));

mu = single(mu(idx, :));
if ~isempty(Sigma),    Sigma    = single(Sigma(idx, idx, :)); end
if ~isempty(invSigma), invSigma = single(invSigma(idx, idx, :)); end
if ~isempty(radius),   radius   = single(radius(idx, :)); end

% Prefer invSigma if user asked and present
if opt.PreferInv && ~isempty(invSigma)
    % rebuild Sigma from invSigma for consistency? Not needed—downstream uses inverse.
    % Keep Sigma empty, pass invSigma via out.SigmaInv for completeness.
    out.SigmaInv = invSigma;
    if isempty(Sigma)
        % provide Sigma too for downstream that expects Sigma
        Sigma = zeros(size(invSigma), 'single');
        for k=1:size(invSigma,3)
            Sigma(:,:,k) = single(pinv(double(invSigma(:,:,k))));
        end
    end
end

% ---- pack output ----
out.mu      = mu;                          % C x K
out.Sigma   = Sigma;                       % C x C x K
out.radius  = radius;                      % C x K or []
out.usedStr = char(used);
out.scaleIn = char(scaleIn);
out.channelOrderIn = char(orderIn);
end

function y = tern(cond, a, b)
if cond, y=a; else, y=b; end
end
