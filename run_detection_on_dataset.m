function R = run_detection_on_dataset(inFolder, outFolder, varargin)
% RUN_DETECTION_ON_DATASET  End-to-end folder runner for FlakeHunter detection.
%
% R = run_detection_on_dataset(inFolder, outFolder, ...
%        'ParamsFile',  ...,     % .json/.mat params for the GMM
%        'FlatfieldFile',[],     % optional flatfield (HxWx3 or 1x1x3)
%        'UsedChannels','RGB',   % channel subset passed to GMM loader
%        'AssumeOrder','auto',   % 'auto'|'RGB'|'BGR' for params
%        'ConfThresh',0.6,       % threshold for confMap -> mask
%        'OpenRadiusPx',1, 'CloseRadiusPx',1,   % morphology
%        'MinArea_um2',0, 'AutoAreaPct',[], 'MinSolidity',0, 'KeepTopK',0,
%        'Magnification',[], 'UmPerPixel',[],   % scale overrides
%        'Chi2Confidence',true,  'UseRadiusGating',true, ...
%        'SaveMask',true, 'SaveStats',true, 'SaveConf',false, ...
%        'Pattern','*.png', 'Recursive',false, 'Verbose',true)
%
% INPUTS
%   inFolder   : folder with input images
%   outFolder  : folder to write outputs (created if missing)
%
% OUTPUT
%   R : struct array (per image) with fields:
%       .imagePath, .maskPath, .statsPath, .ok, .err, .nFlakes, .timeSec
%
% PRECEDENCE RULES (per image)
%   1) Explicit function args (UmPerPixel, Magnification) override everything
%   2) Sidecar JSON via detector.readImageMeta (your version) fields
%   3) Filename fallback via detector.stubMetaFromName (called inside readImageMeta)
%   4) Defaults in detector.umPerPixel lookup (if only magnification known)
%
% NOTE
%   - This runner is I/O glue only; the core math lives in +detector.
%   - No toolboxes required.

tStart = tic;

% -------- parse & defaults --------
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('ParamsFile','',@(s)ischar(s)||isstring(s));
ip.addParameter('FlatfieldFile','',@(s)ischar(s)||isstring(s)||isempty(s));
ip.addParameter('UsedChannels','RGB',@(s)ischar(s)||isstring(s));
ip.addParameter('AssumeOrder','auto',@(s)any(strcmpi(s,{'auto','RGB','BGR'})));

ip.addParameter('ConfThresh',0.6,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('OpenRadiusPx',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('CloseRadiusPx',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('MinArea_um2',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('AutoAreaPct',[],@(x)isempty(x)||isnumeric(x));
ip.addParameter('MinSolidity',0,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('KeepTopK',0,@(x)isnumeric(x)&&isscalar(x));

ip.addParameter('Magnification',[],@(x)ischar(x)||isstring(x)||isnumeric(x)||isempty(x));
ip.addParameter('UmPerPixel',[],@(x)isempty(x)||(isscalar(x)&&x>0));
ip.addParameter('Chi2Confidence',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('UseRadiusGating',true,@(b)islogical(b)||ismember(b,[0 1]));

ip.addParameter('Pattern','*.png',@(s)ischar(s)||isstring(s));
ip.addParameter('Recursive',false,@(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('SaveMask',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('SaveStats',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('SaveConf',false,@(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('Verbose',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.parse(varargin{:});
opt = ip.Results;

inFolder  = char(inFolder);
outFolder = char(outFolder);
assert(isfolder(inFolder), 'Input folder not found: %s', inFolder);
if ~isfolder(outFolder), mkdir(outFolder); end

% Resolve files
files = list_files(inFolder, opt.Pattern, opt.Recursive);
if isempty(files)
    warning('No images found in %s (pattern: %s)', inFolder, opt.Pattern);
    R = struct([]); return
end

% Pre-load flatfield if provided (used as-is for all images)
flatfield = [];
if ~isempty(opt.FlatfieldFile) && exist(opt.FlatfieldFile,'file')==2
    flatfield = local_read_flatfield(opt.FlatfieldFile);
end

% Load/normalize GMM once (global model)
assert(~isempty(opt.ParamsFile), 'ParamsFile is required for the runner.');
rawParams = detector.readGmmParams(opt.ParamsFile);
[params, usedStr] = local_select_channels(rawParams, opt.UsedChannels, opt.AssumeOrder);
paramsStruct = struct('mu', single(params.mu), 'Sigma', single(params.Sigma));
if ~isempty(params.radius)
    paramsStruct.radius = single(params.radius);
end
% Output subfolders
if opt.SaveMask,  outMaskDir  = fullfile(outFolder,'masks');   if ~isfolder(outMaskDir),  mkdir(outMaskDir);  end, end
if opt.SaveConf,  outConfDir  = fullfile(outFolder,'conf');    if ~isfolder(outConfDir),  mkdir(outConfDir);  end, end
if opt.SaveStats, outStatsDir = fullfile(outFolder,'stats');   if ~isfolder(outStatsDir), mkdir(outStatsDir); end

R = repmat(struct('imagePath','','maskPath','','statsPath','','ok',false,'err','','nFlakes',0,'timeSec',0), numel(files),1);

if opt.Verbose
    fprintf('[run] %d images | UsedChannels=%s | ConfThresh=%.3f | Params=%s\n', ...
        numel(files), char(usedStr), opt.ConfThresh, char(opt.ParamsFile));
end

% -------- main loop --------
for i=1:numel(files)
    t0 = tic;
    fp = files{i};
    R(i).imagePath = fp;

    try
        % --- read image
        I = imread(fp);
        if ndims(I)==2, I = repmat(I,1,1,3); end
        if size(I,3)>3, I = I(:,:,1:3); end

        % --- read image meta (your tolerant reader)
         [folder, base, ~] = fileparts(fp);
        M = detector.readImageMeta(folder, base);
        if isempty(fieldnames(M))
            M = detector.stubMetaFromName(base);
        end
        % --- resolve scale precedence
        % Explicit > meta.umPerPixel > meta.magnification via umPerPixel lookup > default 1
        if ~isempty(opt.UmPerPixel)
            umPerPx = double(opt.UmPerPixel);
        else
         umPerPx = local_um_per_pixel(M, opt.Magnification);
            if isempty(umPerPx)
                umPerPx = 1;
            end
        end

        % --- preprocess (flatfield optional, keep intensity type stable)
        procOpts = struct('UsedChannels', usedStr, 'BackgroundCeiling', 241);
        if ~isempty(flatfield)
            procOpts.Flatfield = flatfield;
        end
                proc = detector.preprocess(I, procOpts);

         % --- inference (toolbox-free src/+detector implementation)
        inferOpts = struct('Chi2Conf', logical(opt.Chi2Confidence), ...
                           'ReturnBestComp', true, ...
                           'MinChannelGate', logical(opt.UseRadiusGating));
        infOut = detector.inferGMM(proc, paramsStruct, inferOpts);
        confMap  = infOut.confMap;
        labelMap = infOut.labelMap;
        
        % --- postprocess (threshold, morphology, component filters)
        pp = detector.postprocess(confMap, labelMap, struct( ...
            'ConfThresh',   opt.ConfThresh, ...
            'OpenRadiusPx', opt.OpenRadiusPx, ...
            'CloseRadiusPx',opt.CloseRadiusPx, ...
            'MinArea_um2',  opt.MinArea_um2, ...
            'AutoAreaPct',  opt.AutoAreaPct, ...
            'umPerPixel',   umPerPx, ...
            'MinSolidity',  opt.MinSolidity, ...
            'KeepTopK',     opt.KeepTopK ));

        % --- save outputs
        maskPath = ''; statsPath = ''; confPath = '';
        [~,base,~] = fileparts(fp);

        if opt.SaveMask
            maskPath = fullfile(outMaskDir, [base '_mask.png']);
            imwrite(pp.mask, maskPath);
        end
        if opt.SaveConf
            confPath = fullfile(outConfDir, [base '_conf.tif']);
            imwrite(uint8(max(0,min(255,round(255*pp.mask.*confMap)))), confPath);
        end
        if opt.SaveStats
            statsPath = fullfile(outStatsDir, [base '_stats.csv']);
            writetable(pp.stats, statsPath);
            % also dump a .json sidecar for completeness
            jsonPath = fullfile(outStatsDir, [base '_stats.json']);
            S = table2struct(pp.stats);
           meta = struct('imagePath',fp, 'umPerPixel', umPerPx, 'usedChannels', char(usedStr));
            fid = fopen(jsonPath,'w'); fwrite(fid, jsonencode(J),'char'); fclose(fid);
        end

        % --- finalize record
        R(i).maskPath  = maskPath;
        R(i).statsPath = statsPath;
        R(i).ok        = true;
        R(i).nFlakes   = height(pp.stats);
        R(i).timeSec   = toc(t0);

        if opt.Verbose
            fprintf('[%3d/%3d] %s  flakes=%d  um/px=%.4g  (%.2fs)\n', ...
                i, numel(files), base, R(i).nFlakes, umPerPx, R(i).timeSec);
        end

    catch ME
        R(i).ok = false;
        R(i).err = ME.message;
        if opt.Verbose
            warning('[%3d/%3d] %s  ERROR: %s', i, numel(files), fp, ME.message);
        end
    end
end

if opt.Verbose
    dt = toc(tStart);
    fprintf('[run] done in %.2fs (%.2f imgs/s)\n', dt, numel(files)/dt);
end
end

% ===== helpers =====
function L = list_files(root, pat, recursive)
if recursive
    dd = dir(fullfile(root, '**', pat));
else
    dd = dir(fullfile(root, pat));
end
L = cellfun(@(d,n) fullfile(d,n), {dd.folder}, {dd.name}, 'uni',0);
end

function F = local_read_flatfield(p)
[~,~,e] = fileparts(p);
switch lower(e)
    case '.mat'
        S = load(p); k = fieldnames(S); F = S.(k{1});
    otherwise
        F = imread(p);
end
F = double(F);
if ndims(F)==2, F = repmat(F,1,1,3); end
end
function [out, usedStr] = local_select_channels(rawParams, requested, assumeOrder)
% Convert loaded GMM params (JSON or legacy structs) to the subset/order
% expected by detector.inferGMM (C x K in [0,1], RGB-oriented).

usedStr = upper(char(requested));
usedStr = usedStr(ismember(usedStr,'RGB'));
if isempty(usedStr)
    usedStr = 'RGB';
end

if nargin < 3 || isempty(assumeOrder)
    assumeOrder = 'auto';
end

orderPref = upper(string(assumeOrder));
if orderPref=="AUTO"
    if isfield(rawParams,'channelOrderIn') && ~isempty(rawParams.channelOrderIn)
        orderPref = upper(string(rawParams.channelOrderIn));
    elseif isfield(rawParams,'channelOrder') && ~isempty(rawParams.channelOrder)
        orderPref = upper(string(rawParams.channelOrder));
    elseif isfield(rawParams,'mu_bgr')
        orderPref = "BGR";
    else
        orderPref = "RGB";
    end
elseif ~(orderPref=="RGB" || orderPref=="BGR")
    orderPref = "RGB";
end

% Gather channel-first (C x K) arrays and remember their current order.
muCk = [];
SigmaC = [];
radiusCk = [];
orderStr = '';

if isfield(rawParams,'mu_rgb')
    muCk = double(rawParams.mu_rgb.');
    if isfield(rawParams,'Sigma_rgb'),  SigmaC = double(rawParams.Sigma_rgb); end
    if isfield(rawParams,'radius_rgb'), radiusCk = double(rawParams.radius_rgb.'); end
    orderStr = 'RGB';
    if orderPref=="BGR"
        if isfield(rawParams,'mu_bgr')
            muCk = double(rawParams.mu_bgr.');
            if isfield(rawParams,'Sigma_bgr'),  SigmaC = double(rawParams.Sigma_bgr); end
            if isfield(rawParams,'radius_bgr'), radiusCk = double(rawParams.radius_bgr.'); end
        else
            perm = [3 2 1];
            muCk = muCk(perm,:);
            if ~isempty(SigmaC),  SigmaC = SigmaC(perm,perm,:); end
            if ~isempty(radiusCk), radiusCk = radiusCk(perm,:); end
        end
        orderStr = 'BGR';
    end
else
    if isfield(rawParams,'mu')
        muCk = double(rawParams.mu);
    elseif isfield(rawParams,'mu_bgr')
        muCk = double(rawParams.mu_bgr.');
        orderPref = "BGR";
    else
        error('Params struct missing mu/mu_rgb fields.');
    end
    if orderPref=="BGR"
        if isfield(rawParams,'Sigma_bgr')
            SigmaC = double(rawParams.Sigma_bgr);
        elseif isfield(rawParams,'Sigma')
            SigmaC = double(rawParams.Sigma);
        end
        if isfield(rawParams,'radius_bgr')
            radiusCk = double(rawParams.radius_bgr);
            if size(radiusCk,1) ~= size(muCk,1)
                radiusCk = double(rawParams.radius_bgr.');
            end
        elseif isfield(rawParams,'radius')
            radiusCk = double(rawParams.radius);
        end
    else
        if isfield(rawParams,'Sigma')
            SigmaC = double(rawParams.Sigma);
        elseif isfield(rawParams,'Sigma_rgb')
            SigmaC = double(rawParams.Sigma_rgb);
        end
        if isfield(rawParams,'radius')
            radiusCk = double(rawParams.radius);
        elseif isfield(rawParams,'radius_rgb')
            radiusCk = double(rawParams.radius_rgb);
            if size(radiusCk,1) ~= size(muCk,1)
                radiusCk = double(rawParams.radius_rgb.');
            end
        end
    end

    % Ensure muCk is C x K
    if isfield(rawParams,'K') && ~isempty(rawParams.K)
        expectedK = double(rawParams.K);
    else
        expectedK = size(muCk,2);
    end
    if size(muCk,2) ~= expectedK && size(muCk,1) == expectedK
        muCk = muCk.';
    end
    if size(muCk,1) < size(muCk,2) && size(muCk,2) ~= expectedK
        muCk = muCk.';
    end
    orderStr = local_align_order(orderPref, size(muCk,1));
end

% Reorder to canonical RGB order (keeping only available channels).
canonOrder = local_canonical_from(orderStr);
canonOrder = canonOrder(1:min(numel(canonOrder), size(muCk,1)));
permToCanon = local_order_to_target(orderStr, canonOrder);
muCk = muCk(permToCanon,:);
if ~isempty(SigmaC),  SigmaC = SigmaC(permToCanon, permToCanon, :); end
if ~isempty(radiusCk)
    if size(radiusCk,1) == size(muCk,1)
        radiusCk = radiusCk(permToCanon,:);
    else
        radiusCk = radiusCk(:, permToCanon).';
    end
end

% Convert to components-first for channel selection.
muAll = permute(muCk, [2 1]);              % K x C
K = size(muAll,1);
nChan = size(muAll,2);

if numel(canonOrder) ~= nChan
    error('Canonical channel order (%s) mismatch with parameter count (%d).', canonOrder, nChan);
end

if isempty(SigmaC)
    error('GMM params missing covariance matrices (Sigma).');
end
if size(SigmaC,3) ~= K
    error('Sigma third dimension (%d) did not match K (%d).', size(SigmaC,3), K);
end

% Map requested channels into canonical order.
idx = zeros(1, numel(usedStr));
for ii = 1:numel(usedStr)
    pos = find(canonOrder == usedStr(ii), 1);
    assert(~isempty(pos), 'Requested channel %s missing in params (available: %s).', usedStr(ii), canonOrder);
    idx(ii) = pos;
end

mu = permute(muAll(:, idx), [2 1]);
Sigma = zeros(numel(idx), numel(idx), K, 'like', SigmaC);
for k = 1:K
    Sigma(:,:,k) = SigmaC(idx, idx, k);
end

if ~isempty(radiusCk)
    if size(radiusCk,1) ~= nChan
        radiusCk = radiusCk.';
    end
    radius = radiusCk(idx, :);
else
    radius = [];
end

out = struct('mu', single(mu), 'Sigma', single(Sigma));
if ~isempty(radius)
    out.radius = single(radius);
else
    out.radius = [];
end
usedStr = char(usedStr);
end

function um = local_um_per_pixel(meta, fallbackMagnification)
% Resolve µm/px from meta struct or fallback magnification.
um = [];
magVal = [];
if isstruct(meta) && isfield(meta,'magnification') && ~isempty(meta.magnification)
    magVal = meta.magnification;
elseif ~isempty(fallbackMagnification)
    magVal = fallbackMagnification;
end
if isempty(magVal)
    return;
end
if ischar(magVal) || isstring(magVal)
    magStr = char(string(magVal));
    magVal = str2double(magStr);
    if ~isfinite(magVal)
        tok = regexp(magStr, '(\d+(\.\d+)?)', 'tokens', 'once');
        if ~isempty(tok)
            magVal = str2double(tok{1});
        end
    end
end
if ~isfinite(magVal)
    return;
end
try
    idx = detector.magnificationToIndex(double(magVal));
    um = detector.umPerPixel(idx);
catch
    um = [];
end
end

function orderStr = local_align_order(orderHint, count)
% Normalize an order hint (e.g. 'auto','BGR') to a char vector length COUNT.
base = 'RGB';
hint = upper(char(string(orderHint)));
hint = hint(ismember(hint,'RGB'));
orderStr = '';
seen = false(1,3);
for i = 1:numel(hint)
    ch = hint(i);
    idx = find(base==ch, 1);
    if ~isempty(idx) && ~seen(idx)
        orderStr(end+1) = ch; %#ok<AGROW>
        seen(idx) = true;
    end
end
for i = 1:numel(base)
    if numel(orderStr) >= count, break; end
    if ~seen(i)
        orderStr(end+1) = base(i); %#ok<AGROW>
        seen(i) = true;
    end
end
if isempty(orderStr)
    orderStr = base;
end
if numel(orderStr) < count
    % pad by repeating base sequence if necessary (defensive; typically count<=3)
    while numel(orderStr) < count
        orderStr(end+1) = base( mod(numel(orderStr), numel(base)) + 1 ); %#ok<AGROW>
    end
end
orderStr = orderStr(1:min(numel(orderStr), count));
end

function canon = local_canonical_from(orderStr)
% Derive canonical RGB order limited to channels present in orderStr.
orderStr = upper(char(orderStr));
base = 'RGB';
canon = '';
for i = 1:numel(base)
    ch = base(i);
    if any(orderStr == ch)
        canon(end+1) = ch; %#ok<AGROW>
    end
end
if isempty(canon)
    canon = orderStr;
else
    % append any remaining characters (defensive) preserving appearance order
    for i = 1:numel(orderStr)
        ch = orderStr(i);
        if ~any(canon == ch)
            canon(end+1) = ch; %#ok<AGROW>
        end
    end
end
end

function perm = local_order_to_target(orderStr, target)
% Return indices so that orderStr(perm) == target.
orderStr = upper(char(orderStr));
target = upper(char(target));
perm = zeros(1, numel(target));
for i = 1:numel(target)
    ch = target(i);
    idx = find(orderStr == ch, 1);
    if isempty(idx)
        error('Channel %s missing when reordering (available: %s).', ch, orderStr);
    end
    perm(i) = idx;
end
end
