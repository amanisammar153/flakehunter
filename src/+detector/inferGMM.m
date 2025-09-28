function out = inferGMM(proc, params, opts)
%INFERGMM Pixelwise confidence map using GMM (diagonal or full covariance).
% Inputs:
%   proc   : struct from preprocess() with fields X (N x C in [0,1]), usedIdx, usedStr
%   params : struct with fields:
%       .mu       : C x K (means in [0,1] channel domain matching proc.usedIdx)
%       .Sigma    : C x C x K (covariances)
%       .radius   : C x K (per-channel abs distance gates in [0,1] units) [optional]
%   opts:
%       .Chi2Conf : logical (default true) -> convert d2 to conf via chi2cdf
%       .ReturnBestComp : logical (default true) -> labelMap with best k (0=bg)
%       .MinChannelGate : logical (default true) -> apply per-channel radius gating
%
% Outputs:
%   out.confMap  : HxW single in [0,1]
%   out.labelMap : HxW uint16, 0=background, 1..K=component with min d2 (if gated)
%   out.d2min    : HxW single (min Mahalanobis distance^2)
arguments
    proc struct
    params struct
    opts.Chi2Conf (1,1) logical = true
    opts.ReturnBestComp (1,1) logical = true
    opts.MinChannelGate (1,1) logical = true
end

X = proc.X;            % N x C
C = size(X,2);
N = size(X,1);

mu    = single(params.mu);           % C x K
Sigma = single(params.Sigma);        % C x C x K
if isfield(params,'radius') && ~isempty(params.radius)
    radius = single(params.radius);  % C x K
else
    radius = inf(C, size(mu,2), 'single');
end

% sanity checks
K = size(mu,2);
assert(all(size(mu,1)==C), 'mu dimensions mismatch C.');
assert(all(size(Sigma,1)==C && size(Sigma,2)==C && size(Sigma,3)==K), 'Sigma dims mismatch.');

% precompute inverses and log-dets (stable)
iSigma = zeros(C,C,K,'single');
for k=1*K
    S = double(Sigma(:,:,k));
    % add tiny jitter for stability
    S = S + 1e-9*eye(C);
    iSigma(:,:,k) = single(pinv(S));
end

% per-channel gate mask N x K
if opts.MinChannelGate
    % check absolute channel deviation <= radius for all used channels
    gate = true(N,K);
    for k=1:K
        for c=1:C
            dv = abs( X(:,c) - mu(c,k) );
            gate(:,k) = gate(:,k) & (dv <= radius(c,k));
        end
    end
else
    gate = true(N,K);
end

% compute Mahalanobis d^2 to each component
d2 = inf(N,K,'single');
for k=1:K
    diff = X - mu(:,k).';
    % (x-mu) * invS * (x-mu)'  -> rowwise
    M = iSigma(:,:,k);
    t = sum((diff * M) .* diff, 2); % N x 1
    d2(:,k) = single(t);
end

% best component and confidence
[d2min, kmin] = min(d2, [], 2, 'omitnan');

% gate: set to inf where channel gate fails (guard kmin==0)
bad = false(N,1);
idx = find(kmin>0);                 % only valid component indices
if ~isempty(idx)
    bad(idx) = ~gate(sub2ind(size(gate), idx, kmin(idx)));
end
d2min(bad) = inf;
kmin(bad)  = 0; % background


% convert to confidence [0,1]
if opts.Chi2Conf
    % chi-square with DOF = C
    % conf = 1 - cdf('chi2', d2min, C)
    % Avoid Statistics Toolbox with a robust approx (regularized lower gamma).
    conf = chi2_survival_approx(d2min, C);
else
    % simple inverse distance mapping
    conf = 1 ./ (1 + d2min);
    conf(isinf(d2min)) = 0;
end

% reshape back to HxW
H = size(proc.Iproc,1);
W = size(proc.Iproc,2);
out.confMap  = reshape(single(conf), H, W);
out.d2min    = reshape(single(d2min), H, W);
if opts.ReturnBestComp
    out.labelMap = reshape(uint16(kmin), H, W);
else
    out.labelMap = zeros(H,W,'uint16');
end
out.usedStr = proc.usedStr;

end

function s = chi2_survival_approx(x, v)
% Regularized upper incomplete gamma Q(v/2, x/2), vectorized, toolbox-free.
% Returns ~1 - chi2cdf(x,v).
x = single(x);
v = single(v);
s = zeros(size(x),'single');
finiteMask = isfinite(x);
xf = max(0, double(x(finiteMask)));
vf = double(v);
% Use a continued-fraction for the regularized gamma Q(a, x) with a=v/2.
a = vf/2;
z = xf/2;
% Lentz's algorithm for continued fraction of the incomplete gamma (Q).
% Reference: Numerical Recipes (robust for our v<=3 case typically).
eps1 = 1e-8;
F = zeros(size(z));
for i=1:numel(z)
    zi = z(i);
    if zi==0
        F(i) = 1.0;
        continue
    end
    % Initialization
    b0 = 0;
    C = 1/realmin;
    D = 0;
    f = 1;
    for m=1:200
        a_m = m - a;
        if m==1
            b = 1 + zi - a;
        else
            b = 2*m - 1 + zi - a;
        end
        D = b + a_m * D;
        if D==0, D = realmin; end
        C = b + a_m / C;
        if C==0, C = realmin; end
        D = 1/D;
        delta = C*D;
        f = f * delta;
        if abs(delta-1) < eps1
            break
        end
    end
    % Q(a,z) ≈ exp(-zi + a*log(zi) - gammaln(a)) * f
    F(i) = exp(-zi + a*log(zi) - gammaln(a)) * f;
end
s(finiteMask) = single(max(0,min(1,F)));
% clip infinities to 0
s(~finiteMask) = 0;
end
