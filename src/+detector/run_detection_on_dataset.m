function out = run_detection_on_dataset(imageDir, varargin)
%RUN_DETECTION_ON_DATASET  End-to-end batch detection with per-image magnification.
%   out = run_detection_on_dataset("path/to/tiles", 'Material',"graphene", ...)

tStart = tic;

% -------- parse options
ip = inputParser;
ip.addRequired('imageDir', @(s)ischar(s)||isstring(s));

% params or explicit P
ip.addParameter('Material',"",@(s)ischar(s)||isstring(s));
ip.addParameter('Thickness',"",@(s)ischar(s)||isstring(s));
ip.addParameter('Magnification',[],@(x)isscalar(x)&&isnumeric(x));
ip.addParameter('Params',[],@(s)isstruct(s)||isempty(s));

% detector knobs (defaults can be overridden)
ip.addParameter('UsedChannels',"RGB");
ip.addParameter('StdThresh',4.0);
ip.addParameter('ConfThresh',0.55);
ip.addParameter('SizeThreshUm2',200);

% IO
ip.addParameter('OutDir', fullfile(fh_root(),"results"));
ip.addParameter('SaveOverlays',true,@islogical);
ip.addParameter('SaveMasks',false,@islogical);
ip.addParameter('SaveCrops',true,@islogical);
ip.addParameter('Overwrite',true,@islogical);
ip.addParameter('CleanOutDir',true,@islogical);

% overlays
ip.addParameter('Draw',{'bbox','contour'});
ip.addParameter('LineWidth',2);
ip.addParameter('Alpha',0.25);
ip.addParameter('ShowLabels',true);
ip.addParameter('LabelField','confidence');

% execution
ip.addParameter('Parallel',false,@islogical);
ip.addParameter('Verbose',true,@islogical);

% preprocess/gating extras
ip.addParameter('MaxBackgroundValue',241);
ip.addParameter('RadiusScale',1.0);

% postprocess defaults
ip.addParameter('MorphRadiusPx',1);
ip.addParameter('MinSolidity',0.85);
ip.addParameter('ApplyConfFilter',true);
ip.addParameter('ConfRegionThresh',0.60);

% overview/map (optional)
ip.addParameter('OverviewImage',"",@(s)ischar(s)||isstring(s));
ip.addParameter('OverviewOut',"results/overview_marked.png");
ip.addParameter('XRangeMM',105);
ip.addParameter('YRangeMM',103.333);

ip.addParameter('AutoOverviewIfMissing', true,  @islogical);
ip.addParameter('OverviewSizePx',        [1600 1600], @(v)isnumeric(v)&&numel(v)==2);


% crops
ip.addParameter('CropMarginPx',8);

ip.parse(imageDir, varargin{:});
opt = ip.Results;

% ensure OverviewOut lives under OutDir when not absolute
if ~isfield(opt,'OverviewOut') || isempty(opt.OverviewOut) || strlength(string(opt.OverviewOut))==0
    opt.OverviewOut = fullfile(opt.OutDir, "overview_marked.png");
else
    ov = string(opt.OverviewOut);
    % if it's just a file name (no folder), place it under OutDir
    if isempty(fileparts(char(ov)))
        opt.OverviewOut = fullfile(opt.OutDir, ov);
    end
end

imageDir = char(opt.imageDir);
% Resolve imageDir relative to project root if the given path doesn't exist
if ~isfolder(imageDir)
    maybe = fullfile(fh_root(), string(imageDir));
    if isfolder(maybe)
        imageDir = maybe;
    end
end

assert(isfolder(imageDir), 'Image folder not found: %s', imageDir);

% -------- list images
files = detector.listImages(imageDir);
if isempty(files)
    error('No images found in %s', imageDir);
end
if exist('detector.natsortfiles','file'), files = detector.natsortfiles(files); end

% -------- outputs
outBase = char(opt.OutDir);
jsonDir  = fullfile(outBase,"json");
olDir    = fullfile(outBase,"overlays");
maskDir  = fullfile(outBase,"masks");
cropsDir = fullfile(outBase,"crops");

if opt.CleanOutDir && isfolder(outBase)
    try, rmdir(outBase,'s'); catch, warning('Could not clean %s', outBase); end
end
if ~exist(outBase,'dir'), mkdir(outBase); end
if ~exist(jsonDir,'dir'), mkdir(jsonDir); end
if opt.SaveOverlays && ~exist(olDir,'dir'), mkdir(olDir); end
if opt.SaveMasks    && ~exist(maskDir,'dir'), mkdir(maskDir); end
if opt.SaveCrops    && ~exist(cropsDir,'dir'), mkdir(cropsDir); end

% -------- base params (used as defaults & for 1st mag)
if isempty(opt.Params)
    assert(strlength(opt.Material)>0 && strlength(opt.Thickness)>0 && ~isempty(opt.Magnification), ...
        'Provide Material/Thickness/Magnification or pass Params.');
    Pbase = detector.loadParams(opt.Material, opt.Thickness, opt.Magnification, ...
         'UsedChannels', opt.UsedChannels, ...
         'StdThresh',    opt.StdThresh, ...
         'ConfThresh',   opt.ConfThresh, ...
         'SizeThreshUm2',opt.SizeThreshUm2);
else
    Pbase = opt.Params;
    Pbase.usedChannels  = string(opt.UsedChannels);
    Pbase.stdThresh     = opt.StdThresh;
    Pbase.confThresh    = opt.ConfThresh;
    Pbase.sizeThreshUm2 = opt.SizeThreshUm2;
end

% -------- (optional) parallel
usePar = false;
if opt.Parallel
    try, pc = gcp('nocreate'); if isempty(pc), parpool; end, usePar = true; catch, end
end

% persistent param cache across images (by magnification) — non-nested helper
getPforMag = @(magImg, meta) getPforMag_cached(magImg, meta, Pbase, opt);

% -------- loop
N = numel(files);
summary = repmat(struct('image','', 'n_flakes',0, 'mean_conf',NaN, ...
                        'median_area_um2',NaN, 'json','', 'overlay','', 'mask',''), N, 1);

allFlakesMM = repmat(struct('flake_position_x_mm',0,'flake_position_y_mm',0, ...
                            'id',0,'class_id',0,'confidence',0), 0,1);

if usePar
    tmpSum = cell(N,1); tmpMM = cell(N,1);
    parfor i = 1:N
        [s, mm] = processOne(files{i}, getPforMag, jsonDir, olDir, maskDir, cropsDir, opt); % %% CHANGED
        tmpSum{i} = s; tmpMM{i} = mm(:);
    end
    summary = vertcat(tmpSum{:});
    allFlakesMM = vertcat(tmpMM{:});
else
    for i = 1:N
        if opt.Verbose, fprintf('[%d/%d] %s\n', i, N, files{i}); end
        [summary(i), mm] = processOne(files{i}, getPforMag, jsonDir, olDir, maskDir, cropsDir, opt); % %% CHANGED
        allFlakesMM = [allFlakesMM; mm(:)]; %#ok<AGROW>
    end
end

% -------- pack output
out.summaryTable = struct2table(summary);
out.nImages = N;
out.outDirs = struct('json',jsonDir,'overlays',olDir,'masks',maskDir,'crops',cropsDir);

% optional overview mark-up
% overview mark-up (uses real image if provided, otherwise creates a blank canvas)
% overview mark-up: real image if given, otherwise a blank canvas
makeBlank = true;
if strlength(opt.OverviewImage)>0 && isfile(opt.OverviewImage)
    try
        Iov = imread(opt.OverviewImage);
        makeBlank = false;
    catch
        warning('Failed to read OverviewImage: %s. Using blank canvas.', char(opt.OverviewImage));
    end
end
if makeBlank && opt.AutoOverviewIfMissing
    sz = double(opt.OverviewSizePx);
    Hov = max(256, round(sz(1))); Wov = max(256, round(sz(2)));
    Iov = uint8(255*ones(Hov, Wov, 3));  % white canvas
end

if exist('Iov','var')
    Iov = detector.markOnOverviewMap(Iov, allFlakesMM, ...
            'XRangeMM', opt.XRangeMM, 'YRangeMM', opt.YRangeMM, ...
            'XOffsetMM', getfieldOr(Pbase,'xOffset_mm',0), ...
            'YOffsetMM', getfieldOr(Pbase,'yOffset_mm',0), ...
            'ShowLabels', true, 'LineWidth', 2);
    try
        imwrite(Iov, char(opt.OverviewOut));
    catch ME
        warning('Could not write OverviewOut (%s): %s', char(opt.OverviewOut), ME.message);
    end
else
    warning('Overview not written: no OverviewImage provided and AutoOverviewIfMissing=false.');
end



if opt.Verbose
    fprintf('Done %d images in %.2fs. Total flakes: %d\n', ...
        N, toc(tStart), nansum(out.summaryTable.n_flakes));
end
end

% ==================== worker ====================
function [s, flakes_mm] = processOne(imgPath, getPforMag, jsonDir, olDir, maskDir, cropsDir, opt) % %% CHANGED signature
[folder, base, ext] = fileparts(imgPath);
I = imread(imgPath);

% meta (best-effort)
% meta (best-effort)
meta = detector.readImageMeta(folder, base);

% Phase-1: if no sidecar JSON, infer minimal fields from filename
if isempty(meta) || (isstruct(meta) && isempty(fieldnames(meta)))
    meta = detector.stubMetaFromName(base);   % <-- uses filename patterns
end
if isempty(meta), meta = struct(); end        % ensure struct

% fill core fields
meta.image_base     = string([base ext]);
meta.image_path     = string(fullfile(folder, [base ext]));
meta.image_rel_path = meta.image_base;
if ~isfield(meta,'chip_id')  || isempty(meta.chip_id),  meta.chip_id  = base; end
if ~isfield(meta,'image_id') || isempty(meta.image_id), meta.image_id = base; end


% %% NEW: detect magnification from filename (_xNN)
magImg = parseMagFromName(base);
% if not present in name, try meta; else fallback in getPforMag()

% %% NEW: get per-image params (cached; loads if needed)
Pimg = getPforMag(magImg, meta);

% allow per-run scaling of color-radius
if isfield(Pimg,'gmm') && isfield(Pimg.gmm,'radius_rgb')
    Pimg.gmm.radius_rgb = Pimg.gmm.radius_rgb * opt.RadiusScale;
end

% preprocess (use Pimg flatfield & mean)
Icorr = detector.preprocess(I, Pimg.flatfield, ...
         'FlatfieldMean', Pimg.flatfieldMean, ...
         'MaxBackgroundValue', opt.MaxBackgroundValue);

% infer (adaptive to tame outliers), using Pimg thresholds/channels
[mask, labels, conf] = detector.adaptiveInfer(Icorr, Pimg, 'UseRadiusGating', true);

% postprocess (use Pimg µm/px & size threshold)
[flakes, bwClean, ~] = detector.postprocess(Icorr, mask, labels, conf, Pimg, ...
    'MorphRadiusPx', opt.MorphRadiusPx, ...
    'MinSolidity',   opt.MinSolidity, ...
    'ApplyConfFilter', opt.ApplyConfFilter, ...
    'ConfRegionThresh', opt.ConfRegionThresh);

% export JSON (pass Pimg so JSON reflects correct mag/µm/px)
J = detector.exportFlakes(flakes, meta, Pimg, size(Icorr,[1 2]), jsonDir, 'Pretty', true);
jsonPath = J.jsonPath;

% overlays
overlayPath = "";
if opt.SaveOverlays
    overlayPath = fullfile(olDir, sprintf('%s_overlay.png', base));
    detector.markOnOverview(Icorr, flakes, ...
        'Mask', bwClean, 'Draw', opt.Draw, 'LineWidth', opt.LineWidth, ...
        'Alpha', opt.Alpha, 'ShowLabels', opt.ShowLabels, 'LabelField', opt.LabelField, ...
        'SavePath', overlayPath);
end

% mask
maskPath = "";
if opt.SaveMasks
    maskPath = fullfile(maskDir, sprintf('%s_mask.png', base));
    imwrite(uint8(bwClean)*255, maskPath);
end

% crops
if opt.SaveCrops && ~isempty(flakes)
    detector.saveFlakeCrops(Icorr, flakes, cropsDir, ...
        'MarginPx', opt.CropMarginPx, 'Prefix', base + "_");
end

% absolute stage mm for overview (use Pimg µm/px & offsets)
[H,W,~] = size(Icorr); %#ok<ASGLU>
px2mm = Pimg.umPerPx/1000;
mx0 = getfieldOr(meta,'motor_pos_x_mm',0);
my0 = getfieldOr(meta,'motor_pos_y_mm',0);
xOff = getfieldOr(Pimg,'xOffset_mm',0);
yOff = getfieldOr(Pimg,'yOffset_mm',0);
flakes_mm = repmat(struct('flake_position_x_mm',0,'flake_position_y_mm',0, ...
                          'id',0,'class_id',0,'confidence',0), 0,1);
for k=1:numel(flakes)
    c = flakes(k).centroid_px;
    dx_mm = (c(1) - W/2) * px2mm;
    dy_mm = -(c(2) - H/2) * px2mm; % image y-down, stage y-up
    flakes_mm(end+1,1) = struct( ...
        'flake_position_x_mm', mx0 + dx_mm + xOff, ...
        'flake_position_y_mm', my0 + dy_mm + yOff, ...
        'id', flakes(k).id, 'class_id', flakes(k).class_id, ...
        'confidence', flakes(k).confidence); %#ok<AGROW>
end

% summary row
s = struct();
s.image = imgPath;
s.n_flakes = numel(flakes);
s.mean_conf = meanOrNaN([flakes.confidence]);
s.median_area_um2 = medianOrNaN([flakes.area_um2]);
s.json = jsonPath; s.overlay = overlayPath; s.mask = maskPath;
end

% ==================== helpers ====================
function Pimg = getPforMag_cached(magImg, meta, Pbase, opt)
% Cached parameter loader keyed by magnification (non-nested; parfor-safe)
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType','char','ValueType','any');
end

% decide which magnification to use
if isempty(magImg) || isnan(magImg)
    if isstruct(meta) && isfield(meta,'magnification') && ~isempty(meta.magnification)
        magUse = double(meta.magnification);
    else
        magUse = Pbase.magnification; % fallback to base
    end
else
    magUse = double(magImg);
end

% fast path: reuse base if same mag
if isfield(Pbase,'magnification') && Pbase.magnification == magUse
    Pimg = Pbase;
    return;
end

% build cache key (material|thickness|channels|mag|std|conf|area)
baseKey = sprintf('%s|%s|%s', string(Pbase.material), string(Pbase.chipThickness), string(Pbase.usedChannels));
key = sprintf('%s|%d|std%.3f|conf%.3f|area%.1f', baseKey, magUse, Pbase.stdThresh, Pbase.confThresh, Pbase.sizeThreshUm2);

% fetch or load
if cache.isKey(key)
    Pimg = cache(key);
    return;
end

Pimg = detector.loadParams(Pbase.material, Pbase.chipThickness, magUse, ...
            'UsedChannels', Pbase.usedChannels, ...
            'StdThresh',    Pbase.stdThresh, ...
            'ConfThresh',   Pbase.confThresh, ...
            'SizeThreshUm2',Pbase.sizeThreshUm2);

cache(key) = Pimg;
end



function mag = parseMagFromName(base)
% parse a trailing _xNN or ..._xNN_ pattern
mag = NaN;
tok = regexp(base, '_x(\d+)$', 'tokens', 'once');
if isempty(tok)
    tok = regexp(base, '_x(\d+)_', 'tokens', 'once');
end
if ~isempty(tok)
    mag = str2double(tok{1});
end
end

function v = getfieldOr(S,name,def)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
function v = meanOrNaN(x), if isempty(x), v = NaN; else, v = mean(x); end, end
function v = medianOrNaN(x), if isempty(x), v = NaN; else, v = median(x); end, end
