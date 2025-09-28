function [mask, labelMap, confMap, d2min] = inferGMM(Icorr, G, varargin)
% DETECTOR.INFERGMM  Pixelwise inference from a color GMM.
%
% [mask,labelMap,confMap,d2min] = detector.inferGMM(Icorr, G, ...)
%
% Inputs
%   Icorr : HxWxC image AFTER preprocess (uint8/uint16/single/double), C=1..3
%   G     : struct of GMM params. Supported field names:
%           - mu (C x K)            OR  mu_rgb (K x 3)
%           - Sigma (C x C x K)     OR  Sigma_rgb (3 x 3 x K)
%           - invSigma (C x C x K)  (optional; if absent computed)
%           - radius (C x K)        OR  radius_rgb (K x 3)  (optional)
%           (Values are in the SAME channel space and order as Icorr)
%
% Name-Value options
%   'UsedChannels'      char subset of 'RGB' (default auto from Icorr size)
%   'UseRadiusGating'   logical, default true (per-channel abs gate)
%   'ReturnLabels0ForBG'logical, default true (0 for background)
%   'Chi2Confidence'    logical, default true (use chi2cdf on d2)
%   'ConfThresh'        scalar in [0,1], default 0.5 (used to make mask)
%
% Outputs
%   mask     : logical HxW, confMap >= ConfThresh
%   labelMap : uint16 HxW, best component id (0 if gated or bg)
%   confMap  : single HxW in [0,1]
%   d2min    : single HxW, min Mahalanobis distance^2 (Inf where gated out)

% ---------- parse options ----------
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('UsedChannels','',@(s) ischar(s)||isstring(s));
ip.addParameter('UseRadiusGating',true,@(b) islogical(b)||isscalar(b));
ip.addParameter('ReturnLabels0ForBG',true,@(b) islogical(b)||isscalar(b));
ip.addParameter('Chi2Confidence',true,@(b) islogical(b)||isscalar(b));
ip.addParameter('ConfThresh',0.5,@(x) isnumeric(x)&&isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

% ---------- normalize image to [0,1] single features ----------
I = Icorr;
[H,W,Cin] = size(I);
if isempty(opt.UsedChannels)
    if Cin==1, used = 'R'; else, used = 'RGB'; used = used(1:Cin); end
else
    used = upper(char(opt.UsedChannels));
    used = used(ismember(used,'RGB'));
    if isempty(used), used = 'RGB'; used = used(1:Cin); end
end
C = numel(used);

% Build channel index map R=1,G=2,B=3 relative to I (assumes I is RGB order)
map = struct('R',1,'G',2,'B',3);
idx = arrayfun(@(ch) map.(char(ch)), num2cell(used));
if any(idx > Cin)
    error('UsedChannels=%s but image has only %d channels.', used, Cin);
end

% Extract and scale features to [0,1] (robust to different input classes)
Iu8max = double(intmax(class(uint8(1))));
switch class(I)
    case 'uint8',  X = single(reshape(I(:,:,idx),[],C)) / 255;
    case 'uint16', X = single(reshape(I(:,:,idx),[],C)) / double(intmax('uint16'));
    otherwise,     X = single(reshape( max(0,min(1,double(I(:,:,idx))/255)), [], C ));
end
N = size(X,1);

% ---------- normalize/reshape GMM parameters ----------
[mu, Sigma, invSigma, radius] = normalizeGMMParams(G, C);

K = size(mu,2);
if isempty(invSigma)
    invSigma = zeros(C,C,K,'single');
    for k = 1:K
        S = double(Sigma(:,:,k)) + 1e-9*eye(C);
        % Prefer Cholesky to avoid the pseudoinverse when the covariance is SPD.
    [R,p] = chol(S);
    if p==0
        I = eye(C);
        iSigma(:,:,k) = single(R \ (R' \ I));
    else
        % Fall back to pseudoinverse if the component covariance is ill-conditioned.
        iSigma(:,:,k) = single(pinv(S));
    end
    end
end
if isempty(radius)
    radius = inf(C,K,'single');
end

% ---------- per-channel radius gate ----------
if opt.UseRadiusGating
    gate = true(N,K);
    for k=1:K
        for c=1:C
            gate(:,k) = gate(:,k) & (abs(X(:,c) - mu(c,k)) <= radius(c,k));
        end
    end
else
    gate = true(N,K);
end

% ---------- Mahalanobis distance^2 ----------
d2 = inf(N,K,'single');
for k=1:K
    diff = X - mu(:,k).';
    M = invSigma(:,:,k);
    d2(:,k) = single(sum((diff * M) .* diff, 2));
end

% Best component
[d2min_vec, kmin] = min(d2, [], 2, 'omitnan');

% Apply gate safely (guard kmin==0)
bad = false(N,1);
idxv = find(kmin>0);
if ~isempty(idxv)
    bad(idxv) = ~gate(sub2ind(size(gate), idxv, kmin(idxv)));
end
d2min_vec(bad) = inf;
kmin(bad) = 0;

% ---------- Confidence map ----------
if opt.Chi2Confidence
    conf_vec = chi2_survival(d2min_vec, C); % 1 - chi2cdf(d2, C)
else
    conf_vec = 1 ./ (1 + d2min_vec); conf_vec(isinf(d2min_vec)) = 0;
end

% ---------- reshape outputs ----------
confMap = reshape(single(conf_vec), H, W);
d2min   = reshape(single(d2min_vec), H, W);
if opt.ReturnLabels0ForBG
    labelMap = reshape(uint16(kmin), H, W);
else
    % 1..K everywhere; bg stays at previous best, user can mask via conf
    labelMap = reshape(uint16(max(kmin,1)), H, W);
end
mask = confMap >= single(opt.ConfThresh);

end % inferGMM

% ===== Helpers =====
function [mu, Sigma, invSigma, radius] = normalizeGMMParams(G, C)
% Accepts a variety of field names; returns:
%   mu:     C x K (single)
%   Sigma:  C x C x K (single) or []
%   invSigma: C x C x K (single) or []
%   radius: C x K (single) or []
mu = []; Sigma = []; invSigma = []; radius = [];

% Means
if isfield(G,'mu')
    mu = single(G.mu);  % C x K
elseif isfield(G,'mu_rgb')
    mu = single(G.mu_rgb.'); % (K x 3) -> (3 x K)
end
assert(~isempty(mu), 'GMM: missing mu or mu_rgb.');
K = size(mu,2);
assert(size(mu,1)==C, 'mu has %d channels; expected %d.', size(mu,1), C);

% Covariances/inverses
if isfield(G,'invSigma')
    invSigma = single(G.invSigma);
end
if isempty(invSigma)
    if isfield(G,'Sigma')
        Sigma = single(G.Sigma);
    elseif isfield(G,'Sigma_rgb')
        Sigma = single(G.Sigma_rgb);
    end
    assert(~isempty(Sigma), 'GMM: missing Sigma / Sigma_rgb / invSigma.');
    assert(all(size(Sigma,1)==C & size(Sigma,2)==C & size(Sigma,3)==K), 'Sigma dims mismatch.');
else
    assert(all(size(invSigma,1)==C & size(invSigma,2)==C & size(invSigma,3)==K), 'invSigma dims mismatch.');
end

% Radius (optional)
if isfield(G,'radius')
    radius = single(G.radius); % C x K
elseif isfield(G,'radius_rgb')
    radius = single(G.radius_rgb.'); % (K x 3) -> (3 x K)
end
if ~isempty(radius)
    assert(all(size(radius,1)==C & size(radius,2)==K), 'radius dims mismatch.');
end
end

function conf = chi2_survival(d2, dof)
% 1 - chi2cdf(d2, dof) with toolbox fallback
d2 = double(d2);
if exist('chi2cdf','file')==2
    conf = single(1 - chi2cdf(max(d2,0), dof));
else
    % regularized upper incomplete gamma: Q(dof/2, d2/2)
    a = dof/2; x = max(d2,0)/2;
    conf = single(gammainc(x, a, 'upper'));  % MATLAB base supports this form
end
conf(~isfinite(conf)) = 0;
conf = max(0,min(1,conf));
end
