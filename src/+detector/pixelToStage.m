function out = exportFlakes(flakes, meta, P, imSize, outDir, varargin)
%DETECTOR.EXPORTFLAKES  Save per-image JSON of detected flakes.
%   out = detector.exportFlakes(flakes, meta, P, [H W], outDir, Name,Value,...)
%
% Name-Value options
%   'FileName'     : override base name (default: <image_id> or 'flakes')
%   'IncludeSummary': logical, include per-image summary block (default true)
%   'ImageYAxis'   : 'down'|'up' (default 'down')
%   'StageYAxis'   : 'up'|'down' (default 'up')
%   'Pretty'       : pretty-print JSON (default true)
%
% Output
%   out.jsonPath : written JSON path
%   out.nFlakes  : number written
%   out.summary  : struct with image-level stats (if included)

ip = inputParser;
ip.addParameter('FileName', '', @(s)ischar(s)||isstring(s));
ip.addParameter('IncludeSummary', true, @islogical);
ip.addParameter('ImageYAxis','down',@(s) any(strcmpi(s,{'down','up'})));
ip.addParameter('StageYAxis','up',  @(s) any(strcmpi(s,{'down','up'})));
ip.addParameter('Pretty', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

if ~exist(outDir,'dir'), mkdir(outDir); end

% ----- build list of flake dicts
Jlist = cell(1, numel(flakes));
for i=1:numel(flakes)
    Jlist{i} = detector.flakeStruct(flakes(i), meta, P, imSize, ...
                   'ImageYAxis', opt.ImageYAxis, 'StageYAxis', opt.StageYAxis);
end

% ----- package top-level JSON
payload = struct();
payload.version   = "matlab-port-1";
payload.material  = string(P.material);
payload.thickness = string(P.chipThickness);
payload.magnification = P.magnification;

% summary
if opt.IncludeSummary
    payload.summary = struct( ...
        'n_flakes', numel(flakes), ...
        'mean_conf', meanOrNaN([flakes.confidence]), ...
        'median_area_um2', medianOrNaN([flakes.area_um2]), ...
        'microns_per_pixel', P.umPerPx, ...
        'image_h', imSize(1), 'image_w', imSize(2) );
end

payload.flakes = [Jlist{:}];

% ----- choose file name
base = char(ifelse(~isempty(metaField(meta,'image_id')), meta.image_id, "flakes"));
if ~isempty(opt.FileName), base = char(opt.FileName); end
jsonPath = fullfile(outDir, sprintf('%s_flakes.json', base));

% ----- write JSON
txt = jsonencode(payload);
if opt.Pretty
    txt = utils.jsonencodePretty(payload);
end
fid = fopen(jsonPath,'w'); assert(fid>0, 'Cannot open %s', jsonPath);
fwrite(fid, txt, 'char'); fclose(fid);

out = struct('jsonPath', jsonPath, 'nFlakes', numel(flakes));
if opt.IncludeSummary
    out.summary = payload.summary;
end
end

% ---- helpers
function tf = metaField(meta, name)
tf = isstruct(meta) && isfield(meta, name) && ~isempty(meta.(name));
end

function v = ifelse(cond, a, b)
if cond, v = a; else, v = b; end
end

function m = meanOrNaN(v)
if isempty(v), m = NaN; else, m = mean(v); end
end

function m = medianOrNaN(v)
if isempty(v), m = NaN; else, m = median(v); end
end
