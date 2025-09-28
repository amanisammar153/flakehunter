function out = exportFlakes(flakes, meta, P, imSize, outDir, varargin)
%DETECTOR.EXPORTFLAKES  Save per-image JSON of detected flakes.
%   out = detector.exportFlakes(flakes, meta, P, [H W], outDir, Name,Value,...)
%
% Name-Value:
%   'FileName'        : override base name (default: <image_id> or 'flakes')
%   'IncludeSummary'  : logical, include per-image summary block (default true)
%   'ImageYAxis'      : 'down'|'up' (default 'down')
%   'StageYAxis'      : 'up'|'down' (default 'up')
%   'Pretty'          : pretty-print JSON (default true)

ip = inputParser;
ip.addParameter('FileName','',@(s)ischar(s)||isstring(s));
ip.addParameter('IncludeSummary',true,@islogical);
ip.addParameter('ImageYAxis','down',@(s) any(strcmpi(s,{'down','up'})));
ip.addParameter('StageYAxis','up',@(s) any(strcmpi(s,{'down','up'})));
ip.addParameter('Pretty',true,@islogical);
ip.parse(varargin{:});
opt = ip.Results;

if ~exist(outDir,'dir'), mkdir(outDir); end

% ---- build per-flake dicts
Jlist = cell(1,numel(flakes));
for i=1:numel(flakes)
    Jlist{i} = detector.flakeStruct(flakes(i), meta, P, imSize, ...
        'ImageYAxis', opt.ImageYAxis, 'StageYAxis', opt.StageYAxis);
end
payload_flakes = [Jlist{:}];

% ---- top-level payload (stable schema used elsewhere)
payload = struct();
payload.version            = "matlab-port-1";
payload.material           = string(getfieldOr(P,'material',""));
payload.thickness          = string(getfieldOr(P,'chipThickness',""));
payload.magnification      = getfieldOr(P,'magnification',NaN);
payload.used_channels      = string(getfieldOr(P,'usedChannels',""));
payload.thresholds         = struct( ...
    'std',          getfieldOr(P,'stdThresh',NaN), ...
    'conf',         getfieldOr(P,'confThresh',NaN), ...
    'area_um2_min', getfieldOr(P,'sizeThreshUm2',NaN), ...
    'solidity_min', getfieldOr(P,'minSolidity',[]) );
payload.microns_per_pixel  = getfieldOr(P,'umPerPx',NaN);
payload.stage_offsets_mm   = struct('x',getfieldOr(P,'xOffset_mm',0), ...
                                    'y',getfieldOr(P,'yOffset_mm',0));

% image block
img_base = guessImageBase(meta);
payload.image = struct( ...
    'id',       getfieldOr(meta,'image_id',img_base), ...
    'rel_path', getfieldOr(meta,'image_rel_path',img_base), ...
    'size_px',  [imSize(1), imSize(2)] );

% optional summary
if opt.IncludeSummary
    payload.summary = struct( ...
        'n_flakes',        numel(flakes), ...
        'mean_conf',       meanOrNaN(getNumFieldArray(flakes,'confidence')), ...
        'median_area_um2', medianOrNaN(getNumFieldArray(flakes,'area_um2')) );
end

payload.flakes = payload_flakes;

% ---- choose file name
base = char(getfieldOr(meta,'image_id','flakes'));
if ~isempty(opt.FileName), base = char(opt.FileName); end
jsonPath = fullfile(outDir, sprintf('%s_flakes.json', base));

% ---- write JSON (pretty fallback safe)
try
    if opt.Pretty && exist('utils.jsonencodePretty','file')
        txt = utils.jsonencodePretty(payload);
    else
        txt = jsonencode(payload);
    end
catch
    txt = jsonencode(payload);
end
fid = fopen(jsonPath,'w'); assert(fid>0,'Cannot open %s', jsonPath);
fwrite(fid, txt, 'char'); fclose(fid);

out = struct('jsonPath', jsonPath, 'nFlakes', numel(flakes));
if isfield(payload,'summary'), out.summary = payload.summary; end
end

% ---- helpers
function v = getfieldOr(S,name,def)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
function m = meanOrNaN(v), if isempty(v), m = NaN; else, m = mean(v); end, end
function m = medianOrNaN(v), if isempty(v), m = NaN; else, m = median(v); end, end
function v = getNumFieldArray(S, name)
if ~isstruct(S) || isempty(S) || ~isfield(S,name), v = []; return; end
try, v = double([S.(name)]); catch, v = []; end
end
function b = guessImageBase(meta)
b = getfieldOr(meta,'image_base','');
if ~isempty(b), return; end
p = getfieldOr(meta,'image_path','');
if ~isempty(p), [~,n,e] = fileparts(char(p)); b=[n e]; return; end
b = '';
end
