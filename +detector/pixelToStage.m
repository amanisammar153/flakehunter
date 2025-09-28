function XYmm = pixelToStage(XYpx, imSize, varargin)
% DETECTOR.PIXELTOSTAGE  Convert pixel coords → stage coords (mm).
%
% XYmm = detector.pixelToStage(XYpx, imSize, 'UmPerPixel',0.325, ...)
%
% Inputs
%   XYpx   : N×2 [x y] pixel coords (1-based, pixel centers)
%   imSize : [H W] of the image in pixels
%
% Name–Value (all optional; sensible defaults below)
%   'UmPerPixel'     : scalar µm/px (default 1)
%   'FlipY'          : logical (default true)   % image y-down → stage y-up
%   'Origin'         : 'center'|'topleft'|[x0 y0] in px  (default 'center')
%   'ThetaDeg'       : rotation (ccw, degrees) around origin (default 0)
%   'StageOffsetMM'  : [X0 Y0] mm to add after rotation (default [0 0])
%   'AffineUmPerPx'  : 2×2 optional µm/px matrix (anisotropic + shear). If set, overrides UmPerPixel.
%
% Output
%   XYmm   : N×2 [X Y] in mm (stage coords)
%
% Notes
% - Pixel centers are assumed at integer coordinates (MATLAB convention).
% - Typical usage: FlipY=true, Origin='center', ThetaDeg per camera rotation.

% ---------- parse ----------
ip = inputParser; ip.CaseSensitive = false;
ip.addParameter('UmPerPixel', 1, @(x) isscalar(x) && x>0);
ip.addParameter('FlipY', true, @(b) islogical(b) || ismember(b,[0 1]));
ip.addParameter('Origin', 'center', @(o) ischar(o)||isstring(o)||(isnumeric(o)&&numel(o)==2));
ip.addParameter('ThetaDeg', 0, @(x) isscalar(x));
ip.addParameter('StageOffsetMM', [0 0], @(v) isnumeric(v)&&numel(v)==2);
ip.addParameter('AffineUmPerPx', [], @(A) isempty(A) || (isnumeric(A)&&isequal(size(A),[2 2])));
ip.parse(varargin{:});
opt = ip.Results;

% ---------- normalize inputs ----------
XYpx = double(XYpx);
H = double(imSize(1)); W = double(imSize(2));
if ischar(opt.Origin) || isstring(opt.Origin)
    o = lower(string(opt.Origin));
    switch o
        case "center"
            % center at pixel-center mid-point: ((W+1)/2,(H+1)/2)
            org = [ (W+1)/2 , (H+1)/2 ];
        case "topleft"
            org = [ 1 , 1 ];
        otherwise
            error('Unknown Origin "%s".', o);
    end
else
    org = double(opt.Origin(:)).';  % [x0 y0] in px
end

% Shift to origin-centered pixel space
XY0 = XYpx - org;             % [x y] in px, origin at org
if opt.FlipY
    XY0(:,2) = -XY0(:,2);     % y-down → y-up
end

% Scale: µm/px → mm/px (or use affine)
if ~isempty(opt.AffineUmPerPx)
    % Affine in µm/px (2x2) then to mm
    XYmm = (XY0 * opt.AffineUmPerPx.') / 1000;  % (px)·(µm/px) → µm → mm
else
    s = opt.UmPerPixel / 1000;                  % mm/px
    XYmm = XY0 * s;
end

% Rotate about origin (ccw)
th = deg2rad(opt.ThetaDeg);
if th ~= 0
    R = [cos(th) -sin(th); sin(th) cos(th)];
    XYmm = XYmm * R.';
end

% Add stage offset
XYmm = XYmm + double(opt.StageOffsetMM);
end
