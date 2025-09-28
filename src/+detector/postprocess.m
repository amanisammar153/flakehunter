function out = postprocess(confMap, labelMap, opts)
%POSTPROCESS Threshold + morphology + component filtering, with metrics.
% Inputs:
%   confMap  : HxW single [0,1]
%   labelMap : HxW uint16 (0=bg, 1..K component index from inferGMM)
%   opts:
%     .ConfThresh     : [0..1], default 0.5
%     .MinArea_um2    : scalar (µm^2), default 0
%     .AutoAreaPct    : [] or [low high] percent to derive area floor (by component area px)
%     .umPerPixel     : scalar µm/px (for area conversion)
%     .MinSolidity    : [0..1] (0 disables)
%     .KeepTopK       : integer (0 disables)
%     .CloseRadiusPx  : integer (morph close radius), default 0
%     .OpenRadiusPx   : integer (morph open radius),  default 0
%
% Outputs:
%   out.mask      : logical HxW
%   out.stats     : table with per-flake metrics (area_px, area_um2, centroid, solidity, ecc, meanConf, labelK)
%   out.props     : regionprops struct array (for convenience)
arguments
    confMap {mustBeNumeric}
    labelMap {mustBeNumeric}
    opts.ConfThresh    double = 0.5
    opts.MinArea_um2   double = 0
    opts.AutoAreaPct   double = []
    opts.umPerPixel    double = 1
    opts.MinSolidity   double = 0
    opts.KeepTopK      double = 0
    opts.CloseRadiusPx double = 0
    opts.OpenRadiusPx  double = 0
end

H = size(confMap,1); W = size(confMap,2);
confMap = single(confMap);
labelMap = uint16(labelMap);

% --- threshold ---
bw = confMap >= single(opts.ConfThresh);

% --- morphology (fallbacks without IPT) ---
if opts.CloseRadiusPx > 0
    r = max(1, round(opts.CloseRadiusPx));
    se = strel_disk_(r);
    bw = imdilate_(imerode_(bw, se), se);
end
if opts.OpenRadiusPx > 0
    r = max(1, round(opts.OpenRadiusPx));
    se = strel_disk_(r);
    bw = imerode_(imdilate_(bw, se), se);
end

% --- components & base props ---
CC = bwconncomp(bw);
if CC.NumObjects==0
    out.mask = false(H,W);
    out.stats = table();
    out.props = struct('Centroid',{},'Area',{},'Solidity',{},'Eccentricity',{},'MajorAxisLength',{},'MinorAxisLength',{});
    return;
end

% derive area floor (px) if AutoAreaPct is requested
pxPerUm2 = 1 / (opts.umPerPixel.^2);
areaFloorPx = 0;
if ~isempty(opts.MinArea_um2)
    areaFloorPx = max(areaFloorPx, round(opts.MinArea_um2 * pxPerUm2));
end
if ~isempty(opts.AutoAreaPct)
    pr = sort(cellfun(@numel, CC.PixelIdxList));
    lohi = max(0,min(100,opts.AutoAreaPct));
    if ~isempty(pr)
        A = prctile(double(pr), lohi(1));
        areaFloorPx = max(areaFloorPx, round(A));
    end
end

% filter by area and solidity
maskKeep = true(CC.NumObjects,1);
% compute props incrementally (no IPT dependency)
props = regionprops_basic_(CC, confMap);

for i=1:CC.NumObjects
    area_px = props(i).Area;
    if area_px < areaFloorPx
        maskKeep(i) = false; continue;
    end
    if opts.MinSolidity>0 && props(i).Solidity < opts.MinSolidity
        maskKeep(i) = false; continue;
    end
end

% KeepTopK by mean confidence
if opts.KeepTopK>0
    meanConf = [props.MeanConf]';
    [~, order] = sort(meanConf, 'descend');
    keep_idx = false(CC.NumObjects,1);
    take = min(opts.KeepTopK, CC.NumObjects);
    keep_idx(order(1:take)) = true;
    maskKeep = maskKeep & keep_idx;
end

% build final mask and stats
mask = false(H,W);
keptProps = props(maskKeep);
for i=1:numel(keptProps)
    mask(keptProps(i).PixelIdxList) = true;
end

% optional: attach best component label per region (mode of labelMap)
labelK = zeros(numel(keptProps),1,'uint16');
for i=1:numel(keptProps)
    lab = labelMap(keptProps(i).PixelIdxList);
    lab = lab(lab>0);
    if isempty(lab), labelK(i)=0; else
        labelK(i) = mode(lab);
    end
end

% Pack stats table
area_px = arrayfun(@(p) p.Area, keptProps)';
area_um2 = area_px * (opts.umPerPixel^2);
centroids = vertcat(keptProps.Centroid);
solidity  = arrayfun(@(p) p.Solidity, keptProps)';
ecc       = arrayfun(@(p) p.Eccentricity, keptProps)';
maj       = arrayfun(@(p) p.MajorAxisLength, keptProps)';
minax     = arrayfun(@(p) p.MinorAxisLength, keptProps)';
mconf     = arrayfun(@(p) p.MeanConf, keptProps)';

out.mask = mask;
out.props = keptProps;
out.stats = table(area_px, area_um2, centroids, solidity, ecc, maj, minax, mconf, labelK, ...
    'VariableNames', {'area_px','area_um2','centroid_xy','solidity','eccentricity','major_axis','minor_axis','mean_conf','labelK'});
end

% --------- helpers (toolbox-free) ---------

function se = strel_disk_(r)
% binary disk structuring element
[x,y] = meshgrid(-r:r, -r:r);
se = (x.^2 + y.^2) <= r^2;
end

function B = imdilate_(A, se)
B = conv2(single(A), single(se), 'same') > 0;
end

function B = imerode_(A, se)
% erosion via min filter on support
pad = floor(size(se)/2);
A2 = padarray(A, pad, 'replicate');
B = true(size(A));
idx = find(se);
[sy,sx] = ind2sub(size(se), idx);
for k=1:numel(idx)
    B = B & A2( sy(k):sy(k)+size(A,1)-1, sx(k):sx(k)+size(A,2)-1 );
end
end

function P = regionprops_basic_(CC, confMap)
% Minimal set of props: Area, Centroid, ConvexArea (approx), Solidity, Eccentricity, Axes, MeanConf
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

    % second-moment ellipse
    x0 = x - mean(x); y0 = y - mean(y);
    C = cov(double([x0,y0]));
    if any(isnan(C),'all') || any(isinf(C),'all')
        C = eye(2);
    end
    [V,D] = eig(C);
    ev = sort(diag(D),'descend');
    a = 4*sqrt(max(ev(1),eps));
    b = 4*sqrt(max(ev(2),eps));
    P(i).MajorAxisLength = a;
    P(i).MinorAxisLength = b;
    P(i).Eccentricity = sqrt(max(0,1 - (b/a)^2));

    % crude convex area via alpha shape hull on integer grid (fallback)
    try
        k = convhull(double(x), double(y));
        A = polyarea(double(x(k)), double(y(k)));
        P(i).ConvexArea = A;
    catch
        P(i).ConvexArea = P(i).Area; % fallback
    end
    P(i).Solidity = min(1, P(i).Area / max(P(i).ConvexArea, eps));
end
end
