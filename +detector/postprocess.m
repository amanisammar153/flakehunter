function out = postprocess(confMap, labelMap, varargin)
% DETECTOR.POSTPROCESS  Threshold → morphology → CC filtering → stats.
%
% out = detector.postprocess(confMap, labelMap, ...)
%
% Inputs
%   confMap   : HxW single/double in [0,1]
%   labelMap  : HxW uint16 (0=bg, 1..K=component id)  (optional; pass [] if N/A)
%
% Options (struct or name–value):
%   'ConfThresh'     [0..1] scalar, default 0.5
%   'OpenRadiusPx'   integer >=0, default 0      % toolbox-free morphology
%   'CloseRadiusPx'  integer >=0, default 0
%   'MinArea_um2'    scalar >=0, default 0       % area floor in µm^2
%   'AutoAreaPct'    [] or [p] percentile in [0..100], default []
%   'umPerPixel'     scalar >0, default 1        % µm per pixel
%   'MinSolidity'    [0..1], default 0
%   'KeepTopK'       integer >=0, default 0      % by mean confidence
%
% Outputs
%   out.mask    : logical HxW (final kept pixels)
%   out.stats   : table (area_px, area_um2, centroid_xy, solidity, ecc, axes, mean_conf, labelK)
%   out.props   : struct array (lightweight props + PixelIdxList)
%   out.debug   : struct with intermediate info (optional helpers)

% ---------- parse options ----------
opt = struct( ...
    'ConfThresh', 0.5, ...
    'OpenRadiusPx', 0, ...
    'CloseRadiusPx', 0, ...
    'MinArea_um2', 0, ...
    'AutoAreaPct', [], ...
    'umPerPixel', 1, ...
    'MinSolidity', 0, ...
    'KeepTopK', 0 );
if ~isempty(varargin)
    if isstruct(varargin{1})
        user = varargin{1}; fn = fieldnames(user);
        for k=1:numel(fn), opt.(fn{k}) = user.(fn{k}); end
    else
        if mod(numel(varargin),2)~=0, error('Name-value inputs must come in pairs.'); end
        for k=1:2:numel(varargin), opt.(char(varargin{k})) = varargin{k+1}; end
    end
end

% ---------- normalize inputs ----------
H = size(confMap,1); W = size(confMap,2);
confMap = single(confMap);
if isempty(labelMap), labelMap = zeros(H,W,'uint16'); else, labelMap = uint16(labelMap); end

% ---------- threshold ----------
bw = confMap >= single(opt.ConfThresh);

% ---------- morphology (toolbox-free) ----------
if opt.CloseRadiusPx>0
    se = strel_disk_(max(1,round(opt.CloseRadiusPx)));
    bw = imdilate_(imerode_(bw,se), se);
end
if opt.OpenRadiusPx>0
    se = strel_disk_(max(1,round(opt.OpenRadiusPx)));
    bw = imerode_(imdilate_(bw,se), se);
end

% ---------- components ----------
CC = bwconncomp(bw);
if CC.NumObjects==0
    out.mask=false(H,W); out.stats=table(); out.props=struct('PixelIdxList',{}); out.debug=struct('areaFloorPx',0);
    return
end

% ---------- area floors ----------
pxPerUm2 = 1/(opt.umPerPixel.^2);
areaFloorPx = 0;
if ~isempty(opt.MinArea_um2) && opt.MinArea_um2>0
    areaFloorPx = max(areaFloorPx, round(opt.MinArea_um2 * pxPerUm2));
end
if ~isempty(opt.AutoAreaPct)
    sizes = cellfun(@numel, CC.PixelIdxList);
    p = max(0, min(100, opt.AutoAreaPct(1)));
    A = local_prctile(double(sizes), p);
    areaFloorPx = max(areaFloorPx, round(A));
end

% ---------- props (lightweight, toolbox-free) ----------
props = regionprops_basic_(CC, confMap);

% ---------- filtering ----------
keep = true(CC.NumObjects,1);
for i=1:CC.NumObjects
    if props(i).Area < areaFloorPx, keep(i)=false; continue; end
    if opt.MinSolidity>0 && props(i).Solidity < opt.MinSolidity, keep(i)=false; continue; end
end

% KeepTopK by mean confidence
if opt.KeepTopK>0
    meanConf = [props.MeanConf]';
    [~, order] = sort(meanConf, 'descend');
    allow = false(CC.NumObjects,1);
    allow(order(1:min(opt.KeepTopK,numel(order)))) = true;
    keep = keep & allow;
end

% ---------- build final mask ----------
mask = false(H,W);
keptProps = props(keep);
for i=1:numel(keptProps)
    mask(keptProps(i).PixelIdxList) = true;
end

% ---------- label per region (mode over labelMap >0) ----------
labelK = zeros(numel(keptProps),1,'uint16');
for i=1:numel(keptProps)
    lab = labelMap(keptProps(i).PixelIdxList);
    lab = lab(lab>0);
    if isempty(lab), labelK(i)=0; else, labelK(i)=mode(lab); end
end

% ---------- pack stats ----------
area_px   = arrayfun(@(p) p.Area, keptProps)';
area_um2  = area_px * (opt.umPerPixel^2);
centroids = vertcat(keptProps.Centroid);
solidity  = arrayfun(@(p) p.Solidity, keptProps)';
ecc       = arrayfun(@(p) p.Eccentricity, keptProps)';
maj       = arrayfun(@(p) p.MajorAxisLength, keptProps)';
minax     = arrayfun(@(p) p.MinorAxisLength, keptProps)';
mconf     = arrayfun(@(p) p.MeanConf, keptProps)';
varNames = {'area_px','area_um2','centroid_xy','solidity','eccentricity', ...
            'major_axis','minor_axis','mean_conf','labelK'};

out.mask  = mask;
out.props = keptProps;
out.stats = table(area_px, area_um2, centroids, solidity, ecc, maj, minax, mconf, labelK, ...
  'VariableNames', varNames);
out.debug = struct('areaFloorPx',areaFloorPx,'pxPerUm2',pxPerUm2);
end

% ===== toolbox-free helpers =====
function se = strel_disk_(r)
[x,y] = meshgrid(-r:r,-r:r); se = (x.^2 + y.^2) <= r^2;
end
function B = imdilate_(A,se)
B = conv2(double(A), double(se), 'same') > 0;
end
function B = imerode_(A,se)
sumSE = sum(se(:)); B = conv2(double(A), double(se), 'same') >= (sumSE - eps);
end
function q = local_prctile(x, p)
% tiny percentile (linear interp), toolbox-free
x = sort(x(:)); n = numel(x);
if n==0, q=NaN; return; end
p = max(0,min(100,double(p))); pos = 1 + (n-1)*(p/100);
lo=floor(pos); hi=ceil(pos); w=pos-lo;
lo=max(1,min(n,lo)); hi=max(1,min(n,hi));
q = (1-w)*x(lo) + w*x(hi);
end
function P = regionprops_basic_(CC, confMap)
% Minimal props without IPT: Area, Centroid, Solidity (via hull), Ecc, Axes, MeanConf.
P = repmat(struct('PixelIdxList',[],'Area',0,'Centroid',[NaN NaN], ...
    'ConvexArea',0,'Solidity',1,'Eccentricity',0, ...
    'MajorAxisLength',0,'MinorAxisLength',0,'MeanConf',0), CC.NumObjects,1);
for i=1:CC.NumObjects
    pix = CC.PixelIdxList{i};
    P(i).PixelIdxList = pix;
    P(i).Area = numel(pix);
    [y,x] = ind2sub(CC.ImageSize, pix);
    P(i).Centroid = [mean(x), mean(y)];
    P(i).MeanConf = mean(confMap(pix), 'omitnan');
    % ellipse from second moments
    x0 = double(x) - mean(x); y0 = double(y) - mean(y);
    C = cov([x0,y0]); if any(~isfinite(C), 'all'), C = eye(2); end
    [~,D] = eig(C); ev = sort(diag(D),'descend');
    a = 4*sqrt(max(ev(1),eps)); b = 4*sqrt(max(ev(2),eps));
    P(i).MajorAxisLength = a; P(i).MinorAxisLength = b;
    P(i).Eccentricity = sqrt(max(0,1 - (b/a)^2));
    % convex area via hull (fallback if collinear)
    try
        k = convhull(double(x), double(y));
        A = polyarea(double(x(k)), double(y(k)));
        if isfinite(A) && A>0, P(i).ConvexArea = A; else, P(i).ConvexArea = P(i).Area; end
    catch, P(i).ConvexArea = P(i).Area; end
    P(i).Solidity = min(1, P(i).Area / max(P(i).ConvexArea, eps));
end
end
