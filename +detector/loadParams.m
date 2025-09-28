function P = loadParams(varargin)
% DETECTOR.LOADPARAMS  Resolve GMM, flatfield, channels, and scale for detection.
%
% P = detector.loadParams('ParamsFile', ..., 'FlatfieldFile', ..., 'UsedChannels','RGB', ...
%                         'Magnification', 20, 'UmPerPixel', [], 'AssumeOrder','auto')
%
% Inputs (name–value; all optional, but one of ParamsFile/ParamsStruct is needed)
%   'ParamsFile'      : path to .json/.mat with GMM params
%   'ParamsStruct'    : struct already in memory (same as readGmmParams input)
%   'FlatfieldFile'   : path to flatfield image (HxWx3) or 1x1x3
%   'UsedChannels'    : char subset of 'RGB' (default 'RGB')
%   'Magnification'   : numeric or char (e.g., 20 or '20x')  (default [])
%   'UmPerPixel'      : override µm/px; if empty, computed from Magnification via detector.umPerPixel
%   'AssumeOrder'     : 'auto'|'RGB'|'BGR' for params channel order (default 'auto')
%   'PreferInv'       : prefer invSigma if present (default true)
%
% Output (struct)
%   P.GMM      : struct with fields mu (C×K), Sigma (C×C×K), optional radius (C×K), usedStr
%   P.flatfield: numeric HxWx3 (if provided), else []
%   P.usedStr  : e.g., 'RG'
%   P.umPerPx  : scalar µm/px
%   P.info     : struct with paths, scaleIn, channelOrderIn
%
% Notes
% - Case-insensitive file resolution; tries common image extensions for flatfield.
% - If Magnification empty and UmPerPixel empty, defaults to 1 µm/px (warns).
%

% ---- defaults & parse ----
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('ParamsFile','',@(s)ischar(s)||isstring(s));
ip.addParameter('ParamsStruct',[],@(s)isstruct(s) || isempty(s));
ip.addParameter('FlatfieldFile','',@(s)ischar(s)||isstring(s));
ip.addParameter('UsedChannels','RGB',@(s)ischar(s)||isstring(s));
ip.addParameter('Magnification',[],@(x)ischar(x)||isstring(x)||isnumeric(x)||isempty(x));
ip.addParameter('UmPerPixel',[],@(x)isnumeric(x)||isempty(x));
ip.addParameter('AssumeOrder','auto',@(s)any(strcmpi(s,{'auto','RGB','BGR'})));
ip.addParameter('PreferInv',true,@(b)islogical(b)||isscalar(b));
ip.parse(varargin{:});
opt = ip.Results;

% ---- resolve Params (file or struct) ----
if ~isempty(opt.ParamsStruct)
    G = detector.readGmmParams(opt.ParamsStruct, 'UsedChannels', opt.UsedChannels, ...
        'AssumeOrder', opt.AssumeOrder, 'PreferInv', opt.PreferInv);
elseif ~isempty(opt.ParamsFile)
    paramsFile = local_find_file(opt.ParamsFile);
    G = detector.readGmmParams(paramsFile, 'UsedChannels', opt.UsedChannels, ...
        'AssumeOrder', opt.AssumeOrder, 'PreferInv', opt.PreferInv);
else
    error('Provide ParamsFile or ParamsStruct.');
end

% ---- resolve flatfield (optional) ----
flatfield = [];
if ~isempty(opt.FlatfieldFile)
    ffPath = local_find_flatfield(opt.FlatfieldFile);
    flatfield = local_read_flatfield(ffPath);
end

% ---- resolve µm/px ----
umPerPx = [];
if ~isempty(opt.UmPerPixel)
    umPerPx = double(opt.UmPerPixel);
elseif ~isempty(opt.Magnification)
    umPerPx = detector.umPerPixel(opt.Magnification);  % your existing helper
else
    umPerPx = 1;  % safe default; warn so users don’t forget
    warning('loadParams:NoScale','No Magnification/UmPerPixel provided. Using 1 µm/px.');
end

% ---- pack output ----
P.GMM       = rmfield_if(G, {'SigmaInv'}); % keep fields mu, Sigma, radius, usedStr, scaleIn, channelOrderIn
P.flatfield = flatfield;
P.usedStr   = G.usedStr;
P.umPerPx   = umPerPx;
P.info = struct( ...
    'ParamsFile', tern(~isempty(opt.ParamsFile), char(local_abs(opt.ParamsFile)), ''), ...
    'FlatfieldFile', tern(~isempty(opt.FlatfieldFile), char(local_abs(opt.FlatfieldFile)), ''), ...
    'scaleIn', tern(isfield(G,'scaleIn'), G.scaleIn, ''), ...
    'channelOrderIn', tern(isfield(G,'channelOrderIn'), G.channelOrderIn, ''), ...
    'UsedChannels', G.usedStr );

end

% ===== helpers =====
function p = local_find_file(p0)
% case-insensitive file resolver
p0 = char(p0);
if exist(p0,'file')==2, p = p0; return; end
[d, n, e] = fileparts(p0);
if isempty(e) % try common param extensions
    exts = {'.json','.mat'};
else
    exts = {e};
end
if isempty(d), d = pwd; end
L = dir(fullfile(d, '*'));
names = {L.name};
match = '';
for i=1:numel(names)
    nm = names{i};
    for j=1:numel(exts)
        if strcmpi(nm, [n exts{j}])
            match = fullfile(d, nm); break;
        end
    end
    if ~isempty(match), break; end
end
assert(~isempty(match), 'File not found (case-insensitive search): %s', p0);
p = match;
end

function p = local_find_flatfield(p0)
p0 = char(p0);
if exist(p0,'file')==2, p = p0; return; end
[d,n,e] = fileparts(p0);
if isempty(d), d = pwd; end
exts = tern(isempty(e), {'.png','.tif','.tiff','.jpg','.jpeg','.mat'}, {e});
L = dir(fullfile(d,'*'));
match = '';
for i=1:numel(L)
    nm = L(i).name;
    for j=1:numel(exts)
        if strcmpi(nm, [n exts{j}]), match = fullfile(d,nm); break; end
    end
    if ~isempty(match), break; end
end
assert(~isempty(match), 'Flatfield not found: %s', p0);
p = match;
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
if size(F,3)~=3
    error('Flatfield must be HxWx3 or 1x1x3.');
end
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end

function S = rmfield_if(S, f)
for i=1:numel(f)
    if isfield(S,f{i}), S = rmfield(S,f{i}); end
end
end

function p = local_abs(p)
try, p = string(which(char(p))); catch, p = string(p); end
end
