function P = loadParams(material, chipThickness, magnification, varargin)
%LOADPARAMS  Load detection parameters (GMM, flatfield, magnification, thresholds)
%   P = detector.loadParams("graphene","90nm",20)
%
% Optional name-values:
%   'BaseDir'   : where the Parameters/ folder lives (default: project root)
%   'UsedChannels' : e.g. "BGR" (default "BGR", to match repo)
%   'StdThresh' : Mahalanobis sigma threshold (default 3.5)   % tightened
%   'SizeThreshUm2' : min flake area in µm^2 (default 200)
%   'ConfThresh' : min confidence (default 0.80)               % tightened

% -------- inputs / defaults
ip = inputParser;
ip.addParameter('BaseDir', projectRoot(), @(s)isstring(s)||ischar(s));
ip.addParameter('UsedChannels', "BGR", @(s)isstring(s)||ischar(s));
ip.addParameter('StdThresh', 4.0, @isscalar);     % was 5
ip.addParameter('SizeThreshUm2', 200, @isscalar);
ip.addParameter('ConfThresh', 0.60, @isscalar);   % was 0.5
ip.parse(varargin{:});
opt = ip.Results;

P = struct();
P.material      = string(lower(material));
P.chipThickness = string(chipThickness);
P.magnification = magnification;
P.usedChannels  = string(upper(opt.UsedChannels));
P.stdThresh     = opt.StdThresh;
P.sizeThreshUm2 = opt.SizeThreshUm2;
P.confThresh    = opt.ConfThresh;

% -------- resolve parameter paths
paramDir   = fullfile(opt.BaseDir, "Parameters");
gmmPath    = fullfile(paramDir, "GMM_Parameters", sprintf("%s_%s.json", P.material, P.chipThickness));
% fallback if Python repo uses "Contrasts" folder name
if ~isfile(gmmPath)
    alt = fullfile(opt.BaseDir,"Parameters","Contrasts", ...   % fixed baseDir->opt.BaseDir
                   sprintf('%s_%s.json', lower(material), chipThickness));
    if isfile(alt), gmmPath = alt; end
end

flatSearch = sprintf("%s_%s_%dx.png", P.material, P.chipThickness, P.magnification);
flatPath   = findFileCaseInsensitive(fullfile(paramDir, "Flatfields"), flatSearch);
magPath    = fullfile(paramDir, "Scan_Magnification", sprintf("%dx.json", P.magnification));

assert(isfile(gmmPath),    "GMM JSON not found: %s", gmmPath);
assert(isfile(flatPath),   "Flatfield PNG not found matching: %s", flatSearch);
assert(isfile(magPath),    "Scan magnification JSON not found: %s", magPath);

% -------- load files
P.gmm = detector.readGmmParams(gmmPath);     % includes precomputed inv cov, both BGR & RGB order
P.flatfield = imread(flatPath);
P.flatfield = ensureUint8RGB(P.flatfield);
P.flatfieldMean = squeeze(mean(mean(double(P.flatfield),1),2));  % 1x3, in RGB order

mag = jsondecode(fileread(magPath));
P.viewFieldX_mm = mag.view_field_x;
P.viewFieldY_mm = mag.view_field_y;
P.xOffset_mm    = mag.x_offset;
P.yOffset_mm    = mag.y_offset;

% -------- microns-per-pixel (match Utils/conversion_functions.py)
P.magIndex = detector.magnificationToIndex(P.magnification);
P.umPerPx  = detector.umPerPixel(P.magIndex);  % scalar µm/px

% ---- channel selection (letters must come from 'RGB')
P.usedChannels = upper(string(P.usedChannels));
P.usedChannels = char(P.usedChannels);
P.usedChannels = P.usedChannels(ismember(P.usedChannels,'RGB'));  % keep only R/G/B
if isempty(P.usedChannels)
    P.usedChannels = 'RGB';
end

% Optional postprocess defaults (used in Step 3; safe to set now)
if ~isfield(P,'minSolidity')   || isempty(P.minSolidity),    P.minSolidity = 0.8; end
if ~isfield(P,'morphRadiusPx') || isempty(P.morphRadiusPx),  P.morphRadiusPx = 3; end

% -------- channel order bookkeeping
% JSON stores B,G,R; MATLAB images are R,G,B.
P.channelMap.BGR_to_RGB = [3 2 1];   % index mapping to reorder means/covs if needed
P.channelMap.RGB_to_BGR = [3 2 1];

% -------- sanity
assert(all(size(P.flatfield,3)==3), 'Flatfield must be RGB image.');

% --- QA checks (place at the END, after P is fully populated)
assert(isfield(P,'flatfield') && ~isempty(P.flatfield) && size(P.flatfield,3)==3, ...
    'Flatfield missing or not RGB.');
assert(isfield(P,'gmm') && isfield(P.gmm,'K') && P.gmm.K>=1, ...
    'GMM params missing or empty.');
assert(isfield(P,'umPerPx') && P.umPerPx>0, ...
    'Magnification µm/px not set.');

% ensure flatfieldMean exists (keep as 1x3 row)
if ~isfield(P,'flatfieldMean') || isempty(P.flatfieldMean)
    P.flatfieldMean = reshape(squeeze(mean(mean(double(P.flatfield),1),2)), 1, 3);
end

end

% ----- helpers (unchanged)
function r = projectRoot()
here = fileparts(mfilename('fullpath'));
cand = { fullfile(here,".."), fullfile(here,"..","..") };
for i=1:numel(cand)
    p = fullfile(cand{i});
    if isfolder(fullfile(p,"Parameters"))
        r = string(fullfile(p)); return;
    end
end
r = string(pwd);
end

function I = ensureUint8RGB(I)
if ~isa(I,'uint8'), I = im2uint8(I); end
if size(I,3)~=3, error('Expected 3-channel flatfield image.'); end
end

function p = findFileCaseInsensitive(dirPath, fileName)
d = dir(dirPath);
cands = {d.name};
mask = strcmpi(cands, fileName);
idx  = find(mask,1);
if isempty(idx)
    error('File not found (case-insensitive search): %s in %s', fileName, dirPath);
end
p = fullfile(dirPath, cands{idx});
end
