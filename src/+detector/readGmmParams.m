function G = readGmmParams(jsonPath)
%READGMMPARAMS  Load GMM parameters (JSON with components "1","2",...)
% Output fields:
%   .K
%   .mu_bgr      [K x 3]   (B,G,R order from Python)
%   .Sigma_bgr   [3 x 3 x K]
%   .radius_bgr  [K x 3]
%   .mu_rgb, .Sigma_rgb, .radius_rgb  (reordered to MATLAB RGB)
%   .invSigma_*  precomputed inverses (stable)

raw = jsondecode(fileread(jsonPath));

% --- keys may be unordered: enforce numeric order 1..K
keys = fieldnames(raw);
idxNum = nan(numel(keys),1);
for i = 1:numel(keys), idxNum(i) = str2double(keys{i}); end
[~, order] = sort(idxNum, 'ascend');
keys = keys(order);
K = numel(keys);

mu_bgr     = zeros(K,3);
Sigma_bgr  = zeros(3,3,K);
radius_bgr = zeros(K,3);

for k = 1:K
    comp = raw.(keys{k});
    % JSON fields (Python): contrast.{b,g,r}, covariance_matrix (row-major 3x3),
    % color_radius.{b,g,r}. Cast to double to be safe.
    mu_bgr(k,:)      = double([comp.contrast.b, comp.contrast.g, comp.contrast.r]);

    % reshape row-major -> MATLAB (column-major) needs a transpose
    S = reshape(double(comp.covariance_matrix), [3,3]).';
    % Symmetrize to kill tiny asymmetry
    S = 0.5*(S + S.');
    Sigma_bgr(:,:,k) = S;

    if isfield(comp,'color_radius')
        radius_bgr(k,:) = double([comp.color_radius.b, comp.color_radius.g, comp.color_radius.r]);
    else
        radius_bgr(k,:) = [inf inf inf]; % if absent, effectively disables gating
    end
end

% Reorder for MATLAB RGB images
bgr2rgb = [3 2 1];
mu_rgb      = mu_bgr(:, bgr2rgb);
radius_rgb  = radius_bgr(:, bgr2rgb);
Sigma_rgb   = zeros(3,3,K);
for k = 1:K
    Sigma_rgb(:,:,k) = Sigma_bgr(bgr2rgb, bgr2rgb, k);
end

% Stable inverse (prefer Cholesky; fall back to tweak or pinv)
invSigma_bgr = zeros(3,3,K);
invSigma_rgb = zeros(3,3,K);
for k = 1:K
    invSigma_bgr(:,:,k) = invSPD(Sigma_bgr(:,:,k));
    invSigma_rgb(:,:,k) = invSPD(Sigma_rgb(:,:,k));
end

G = struct( ...
    'K', K, ...
    'mu_bgr', mu_bgr, ...
    'Sigma_bgr', Sigma_bgr, ...
    'radius_bgr', radius_bgr, ...
    'mu_rgb', mu_rgb, ...
    'Sigma_rgb', Sigma_rgb, ...
    'radius_rgb', radius_rgb, ...
    'invSigma_bgr', invSigma_bgr, ...
    'invSigma_rgb', invSigma_rgb);

% Final sanity
assert(isfield(G,'mu_rgb')       && isequal(size(G.mu_rgb), [G.K,3]),        'mu_rgb missing/bad size');
assert(isfield(G,'invSigma_rgb') && isequal(size(G.invSigma_rgb),[3,3,G.K]), 'invSigma_rgb missing/bad size');

end

% ---------- helpers ----------
function iS = invSPD(S)
% try Cholesky; if not SPD, nudge diagonal; last resort pinv
S = 0.5*(S + S.');  % symmetrize
epsBase = 1e-9;
for j = 1:3
    [R,p] = chol(S);
    if p==0
        % Solve inv via Cholesky: S = R'R => invS = R \ (R' \ I)
        I = eye(size(S));
        iS = R \ (R' \ I);
        return;
    end
    S = S + epsBase*eye(3);      % nudge and retry
    epsBase = epsBase*10;
end
% fallback
iS = pinv(S);
end
