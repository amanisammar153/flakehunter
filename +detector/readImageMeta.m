function M = readImageMeta(imgPath, varargin)
% DETECTOR.READIMAGEMETA  Tolerant metadata reader for an image.
%
% M = detector.readImageMeta('...\image_123_x20_y-1.23_0.45mm.png', ...)
%
% Looks for a sidecar JSON next to the image (same basename, .json).
% If missing or partial, falls back to parsing the filename via stubMetaFromName.
%
% Name–Value (all optional):
%   'DefaultMag'   : [] or numeric/char (e.g., 20 or '20x')
%   'DefaultUmPerPx' : [] or numeric µm/px
%   'PreferUmPerPx'  : logical, prefer direct µm/px if available (default true)
%
% Output fields (always present; may be empty [] when unknown):
%   .imagePath        : char
%   .imageName        : char
%   .magnification    : numeric [] if unknown
%   .umPerPixel       : numeric [] if unknown
%   .stageXY_mm       : [X Y] mm [] if unknown
%   .usedChannels     : 'RGB' by default (can be overridden from json)
%   .notes            : freeform string (if provided)
%
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('DefaultMag', [], @(x) ischar(x) || isstring(x) || isnumeric(x) || isempty(x));
ip.addParameter('DefaultUmPerPx', [], @(x) isempty(x) || (isscalar(x) && x>0));
ip.addParameter('PreferUmPerPx', true, @(b)islogical(b)||ismember(b,[0 1]));
ip.parse(varargin{:});
opt = ip.Results;

assert(exist(imgPath,'file')==2, 'Image not found: %s', imgPath);
imgPath = char(imgPath);
[folder, base, ~] = fileparts(imgPath);
jsonPath = fullfile(folder, [base '.json']);

M = struct('imagePath',imgPath,'imageName',[base], ...
           'magnification',[],'umPerPixel',[],'stageXY_mm',[], ...
           'usedChannels','RGB','notes','');

% --- 1) read sidecar JSON if present ---
S = struct();
if exist(jsonPath,'file')==2
    try
        raw = fileread(jsonPath);
        S = jsondecode(raw);
    catch ME
        warning('readImageMeta:BadJSON','Failed to read %s (%s). Falling back to filename parsing.', jsonPath, ME.message);
    end
end

% Map common JSON keys (tolerant to naming)
j = lowerFields(S);

% magnification: prefer explicit numeric; support strings like '20x'
mag = [];
if isfield(j,'magnification'), mag = j.magnification; end
if isempty(mag)
    if isfield(j,'mag'), mag = j.mag; end
end
if ischarOrString(mag)
    [magParsed,~] = detector.magnificationToIndex(mag);
    mag = magParsed;
end
if ~isempty(mag) && (~isnumeric(mag) || ~isscalar(mag) || ~isfinite(mag)), mag = []; end

% µm/px
upp = [];
klist = {'umperpixel','um_per_pixel','um_per_px','um_perpix','pixel_size_um'};
for k=1:numel(klist)
    if isfield(j,klist{k}), upp = j.(klist{k}); break; end
end
if ~isempty(upp) && (~isnumeric(upp) || ~isscalar(upp) || ~(upp>0)), upp = []; end

% stage XY (mm)
xy = [];
klist = {'stage','stage_xy','stage_mm','xy_mm'};
for k=1:numel(klist)
    if isfield(j,klist{k}), xy = j.(klist{k}); break; end
end
if ~isempty(xy)
    xy = double(xy(:)).';
    if numel(xy)~=2 || any(~isfinite(xy)), xy = []; end
end

% used channels?
uch = 'RGB';
if isfield(j,'usedchannels'), uch = upper(string(j.usedchannels)); end
if isfield(j,'channels'),      uch = upper(string(j.channels)); end
uch = char(uch);
uch = uch(ismember(uch,'RGB'));
if isempty(uch), uch = 'RGB'; end

% notes
nts = '';
if isfield(j,'notes'), nts = string(j.notes); end

% --- 2) fallback to filename parsing for any missing fields ---
need.mag = isempty(mag);
need.upp = isempty(upp);
need.xy  = isempty(xy);
stub = detector.stubMetaFromName(imgPath);

if need.mag && ~isempty(stub.magnification), mag = stub.magnification; end
if need.upp && ~isempty(stub.umPerPixel),    upp = stub.umPerPixel; end
if need.xy  && ~isempty(stub.stageXY_mm),    xy  = stub.stageXY_mm; end

% --- 3) defaults if still missing ---
if isempty(mag) && ~isempty(opt.DefaultMag)
    [m,~] = detector.magnificationToIndex(opt.DefaultMag);
    mag = m;
end
if isempty(upp) && ~isempty(opt.DefaultUmPerPx)
    upp = double(opt.DefaultUmPerPx);
end
% If still no um/px but we have magnification, compute via lookup
if isempty(upp) && ~isempty(mag)
    try
        upp = detector.umPerPixel(mag);
    catch
        % leave empty
    end
end

% finalize
M.magnification = mag;
M.umPerPixel    = upp;
M.stageXY_mm    = xy;
M.usedChannels  = uch;
M.notes         = char(nts);
end

% ===== helpers =====
function J = lowerFields(S)
J = struct();
if ~isstruct(S) || isempty(fieldnames(S)), return; end
f = fieldnames(S);
for i=1:numel(f)
    J.(lower(f{i})) = S.(f{i});
end
end

function tf = ischarOrString(x)
tf = (ischar(x) || isstring(x));
end
