function M = stubMetaFromName(p)
% DETECTOR.STUBMETAFROMNAME  Heuristic parse of magnification / stage from filename.
%
% Patterns supported (case-insensitive):
%   *_x20*, *20x*, *mag20*, etc.         → magnification = 20
%   *_X=12.34mm_Y=-5.6mm* or *_x12.34_y-5.6mm* → stageXY_mm = [12.34 -5.6]
%   *_upp0.325um* or *_umpp0.325*        → umPerPixel ≈ 0.325
%
% Returns fields (empty if not present):
%   .magnification (numeric)
%   .umPerPixel    (numeric, µm/px)
%   .stageXY_mm    ([X Y] mm)
%
M = struct('magnification',[],'umPerPixel',[],'stageXY_mm',[]);
[~,name,~] = fileparts(char(p));
s = lower(strrep(name,' ',''));

% magnification
patMag = {'(\d+)\s*x','x\s*(\d+)','mag(\d+)','(\d+)mag'};
for i=1:numel(patMag)
    tok = regexp(s, patMag{i}, 'tokens', 'once');
    if ~isempty(tok), M.magnification = str2double(tok{1}); break; end
end

% um/px
tok = regexp(s, 'u(m|mic)pp([0-9]*\.?[0-9]+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(s, 'upp([0-9]*\.?[0-9]+)u?m?', 'tokens', 'once');
end
if ~isempty(tok)
    M.umPerPixel = str2double(tok{end});
end

% stage XY in mm (several variants)
tok = regexp(s, 'x=([-+]?[0-9]*\.?[0-9]+)mm[_-]?y=([-+]?[0-9]*\.?[0-9]+)mm', 'tokens', 'once');
if isempty(tok)
    tok = regexp(s, 'x([-+]?[0-9]*\.?[0-9]+)mm[_-]?y([-+]?[0-9]*\.?[0-9]+)mm', 'tokens', 'once');
end
if isempty(tok)
    tok = regexp(s, 'xy([-+]?[0-9]*\.?[0-9]+),([-+]?[0-9]*\.?[0-9]+)mm', 'tokens', 'once');
end
if ~isempty(tok)
    M.stageXY_mm = [str2double(tok{1}), str2double(tok{2})];
end
end
