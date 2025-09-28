function um = umPerPixel(mag, varargin)
% DETECTOR.UMPERPIXEL  Resolve micrometers-per-pixel (µm/px) for an image.
%
% um = detector.umPerPixel(mag, 'UmPerPixelOverride',[], 'SensorPixelUm',6.5, ...
%                          'ObjectiveMag',[], 'TubeLensMag',1, 'CameraBinning',1, ...
%                          'Lookup',[])
%
% Inputs
%   mag  : numeric or char ('20','20x','x20', etc.). If empty, uses optics or Lookup.
%
% Name–Value options (all optional):
%   'UmPerPixelOverride' : if provided (scalar >0), returned directly.
%   'SensorPixelUm'      : camera pixel size in µm (default 6.5 µm).
%   'ObjectiveMag'       : objective magnification (e.g., 20).
%   'TubeLensMag'        : tube/relay lens total mag factor (default 1).
%   'CameraBinning'      : binning factor (1=none, 2=2×2, etc.) (default 1).
%   'Lookup'             : map/struct for mag→µm/px. If empty, a default table is used.
%
% Output
%   um : scalar µm/px
%
% Precedence:
%   1) UmPerPixelOverride
%   2) Optics formula (if enough info)
%   3) Lookup table using parsed magnification
%
% Notes:
% - Default lookup assumes 6.5 µm pixels, 1× tube lens, 1× binning:
%     5×→1.30, 10×→0.65, 20×→0.325, 50×→0.13
% - You can pass a custom Lookup (struct or containers.Map), e.g. struct('20',0.312,'50',0.125)

% ---------- parse ----------
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('UmPerPixelOverride',[],@(x) isempty(x) || (isscalar(x)&&x>0));
ip.addParameter('SensorPixelUm',6.5,@(x) isscalar(x) && x>0);
ip.addParameter('ObjectiveMag',[],@(x) isempty(x) || (isscalar(x)&&x>0));
ip.addParameter('TubeLensMag',1,@(x) isscalar(x) && x>0);
ip.addParameter('CameraBinning',1,@(x) isscalar(x) && x>0);
ip.addParameter('Lookup',[],@(x) isempty(x) || isstruct(x) || isa(x,'containers.Map'));
ip.parse(varargin{:});
opt = ip.Results;

% 1) direct override
if ~isempty(opt.UmPerPixelOverride)
    um = double(opt.UmPerPixelOverride);
    return
end

% normalise magnification input
magParsed = local_parse_mag(mag);

% 2) optics-based compute (preferred when we know the optics)
objMag = opt.ObjectiveMag;
if isempty(objMag) && ~isempty(magParsed)
    objMag = magParsed;
end

if ~isempty(objMag)
    sysMag = objMag * opt.TubeLensMag;
    um = (opt.SensorPixelUm * opt.CameraBinning) / sysMag;
    return
end

% 3) lookup table (custom or default)
L = opt.Lookup;
if isempty(L)
    % default table: 6.5 µm pixel, 1× tube lens, 1× binning
    L = struct('5',1.30, '10',0.65, '20',0.325, '50',0.13, '100',0.065);
end

if isa(L,'containers.Map')
    key = num2str(magParsed);
    if isKey(L, key), um = L(key); return; end
    error('umPerPixel:LookupMiss','Magnification "%s" not in Lookup map.', key);
else
    % struct lookup (case-insensitive on fieldnames)
    f = fieldnames(L);
    key = num2str(magParsed);
    idx = find(strcmpi(f, key), 1);
    if ~isempty(idx)
        um = L.(f{idx});
        return
    end
    error('umPerPixel:LookupMiss','Magnification "%s" not in Lookup struct.', key);
end
end

% ===== helpers =====
function m = local_parse_mag(mag)
% returns [] if cannot parse
if isempty(mag), m = []; return; end
if isnumeric(mag), m = double(mag); return; end
s = strtrim(string(mag));
% strip any 'x' or 'X' anywhere
s = erase(lower(s), 'x');
% drop spaces
s = regexprep(s,'\s+','');
v = str2double(s);
if isnan(v) || v<=0, m = []; else, m = v; end
end
