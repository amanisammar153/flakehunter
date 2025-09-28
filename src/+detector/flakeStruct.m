function J = flakeStruct(flake, meta, P, imSize, varargin)
%DETECTOR.FLAKESTRUCT  Build one flake's JSON-ready struct.
%   J = detector.flakeStruct(flake, meta, P, [H W], 'ImageYAxis','down','StageYAxis','up')

% ---- options
ip = inputParser;
ip.addParameter('ImageYAxis','down',@(s) any(strcmpi(s,{'down','up'})));
ip.addParameter('StageYAxis','up',  @(s) any(strcmpi(s,{'down','up'})));
ip.parse(varargin{:});
opt = ip.Results;

H = imSize(1); W = imSize(2);
umPerPx = getfieldOr(P,'umPerPx',0);  px2mm = umPerPx/1000;

% identity
chip_id  = char(string(getfieldOr(meta,'chip_id','unknown_chip')));
image_id = char(string(getfieldOr(meta,'image_id','unknown_image')));

% centroid px -> delta mm from image center (y-down image coords)
c   = flake.centroid_px;   % [x y]
dx  = (c(1) - W/2) * px2mm;
dyi = (c(2) - H/2) * px2mm;                 % image downwards
% map image y to stage y according to conventions
imgYsign   = iff(strcmpi(opt.ImageYAxis,'down'), +1, -1);
stageYsign = iff(strcmpi(opt.StageYAxis,'up'),   +1, -1);
dy  = - dyi * imgYsign * stageYsign;

% absolute stage (mm)
mx0 = getfieldOr(meta,'motor_pos_x_mm',0);
my0 = getfieldOr(meta,'motor_pos_y_mm',0);
xOff = getfieldOr(P,'xOffset_mm',0);
yOff = getfieldOr(P,'yOffset_mm',0);
x_mm = mx0 + dx + xOff;
y_mm = my0 + dy + yOff;

% bbox ints
bb = flake.bbox_px; %#ok<NASGU>
bbox = struct('x', round(flake.bbox_px(1)), 'y', round(flake.bbox_px(2)), ...
              'w', round(flake.bbox_px(3)), 'h', round(flake.bbox_px(4)));

% package
J = struct();
J.chip_id    = chip_id;
J.image_id   = image_id;
J.flake_id   = flake.id;
J.class_id   = flake.class_id;
J.confidence = roundSafe(getfieldOr(flake,'confidence',NaN), 6);
J.false_positive_probability = roundSafe(getfieldOr(flake,'false_positive_prob',NaN), 6);

J.bbox         = bbox;
J.area_px      = getfieldOr(flake,'area_px',NaN);
J.area_um2     = roundSafe(getfieldOr(flake,'area_um2',NaN), 3);
J.solidity     = roundSafe(getfieldOr(flake,'solidity',NaN), 6);
J.eccentricity = roundSafe(getfieldOr(flake,'eccentricity',NaN), 6);
J.perimeter_px = roundSafe(getfieldOr(flake,'perimeter_px',NaN), 3);
J.orientation_deg = roundSafe(getfieldOr(flake,'orientation_deg',NaN), 3);
J.major_axis_px   = roundSafe(getfieldOr(flake,'major_axis_px',NaN), 3);
J.minor_axis_px   = roundSafe(getfieldOr(flake,'minor_axis_px',NaN), 3);

% mean colors
mrgb = getvecOr(flake,'mean_rgb',[NaN NaN NaN]);
mbgr = getvecOr(flake,'mean_bgr',[NaN NaN NaN]);
J.mean_contrast_b = roundSafe(mbgr(1),3);
J.mean_contrast_g = roundSafe(mbgr(2),3);
J.mean_contrast_r = roundSafe(mbgr(3),3);
J.mean_rgb        = roundVec(mrgb,3);

J.used_channels       = char(string(getfieldOr(flake,'used_channels','')));
J.microns_per_pixel   = roundSafe(umPerPx, 6);
J.image_rel_path      = char(string(getfieldOr(meta,'image_rel_path','')));

J.motor_pos_x_mm      = roundSafe(mx0, 6);
J.motor_pos_y_mm      = roundSafe(my0, 6);
J.flake_position_x_mm = roundSafe(x_mm, 6);
J.flake_position_y_mm = roundSafe(y_mm, 6);
end

% ---- helpers
function v = getfieldOr(S,name,def)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
function s = iff(cond, a, b), if cond, s=a; else, s=b; end, end
function x = roundSafe(x,n), if isnan(x), return; end, x = round(x,n); end
function v = getvecOr(S,name,def)
v = def;
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    vv = double(S.(name)); if isvector(vv) && numel(vv)>=3, v = vv(:).'; end
end
end
function v = roundVec(v,n), if any(isnan(v)), return; end, v = round(v,n); end
