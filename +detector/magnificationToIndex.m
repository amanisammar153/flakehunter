function [magNum, idx] = magnificationToIndex(mag, presets)
% DETECTOR.MAGNIFICATIONTOINDEX  Parse magnification and map to preset index.
%
% [magNum, idx] = detector.magnificationToIndex(mag, presets)
%
% Inputs
%   mag      : numeric or char (e.g., 20, '20', '20x', 'x20', '  50 X ')
%   presets  : optional numeric vector of allowed mags (default [5 10 20 50 100])
%
% Outputs
%   magNum   : numeric magnification (empty if cannot parse)
%   idx      : index into presets matching magNum (0 if not found or presets omitted)
%
% Notes
% - Matching is exact on value (after parsing). You can extend 'presets' to your set.

if nargin < 2 || isempty(presets)
    presets = [5 10 20 50 100];
else
    validateattributes(presets, {'numeric'}, {'vector','finite','positive'});
end

magNum = local_parse_mag(mag);
if isempty(magNum)
    idx = 0;
    return
end

% exact match in presets (treat them as integers)
tol = 1e-9;
eqv = abs(presets - magNum) < tol;
if any(eqv)
    idx = find(eqv,1,'first');
else
    idx = 0;
end
end

% ---- helpers ----
function v = local_parse_mag(m)
% returns [] if cannot parse a positive number
if isempty(m)
    v = [];
    return
end
if isnumeric(m)
    v = double(m);
    if ~isfinite(v) || v <= 0, v = []; end
    return
end
s = string(m);
s = lower(strtrim(s));
s = erase(s, 'x');            % remove any x/X
s = regexprep(s,'\s+','');    % remove spaces
v = str2double(s);
if isnan(v) || v <= 0
    v = [];
end
end
