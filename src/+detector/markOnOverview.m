function [Iout, meta] = markOnOverview(Irgb, flakes, varargin)
% DETECTOR.MARKONOVERVIEW  Draw flake annotations (bbox/contour/rotated) on an image.
%   [Iout, meta] = detector.markOnOverview(Irgb, flakes, 'Mask', bw, 'Draw', {'bbox','contour'}, ...)
%
% Name-Value (all optional):
%   'Mask'        : HxW logical mask used to draw contours (default: [])
%   'Draw'        : 'bbox' | 'contour' | 'rotated' | cellstr of any (default: {'bbox'})
%   'LineWidth'   : positive scalar (default: 2)
%   'Alpha'       : 0..1 overlay alpha for filled shapes (default: 0.25) [used by fallbacks lightly]
%   'ShowLabels'  : logical (default: true)
%   'LabelField'  : char/string field in flake to display, e.g. 'confidence' (default: 'confidence')
%   'Colormap'    : 'class' or 'fixed' (default: 'class')
%   'FixedColor'  : [R G B] 0..255 used when Colormap='fixed' (default: [0 255 255])
%   'SavePath'    : path to write PNG (default: '')
%
% Works with/without Computer Vision Toolbox. Uses insertShape/insertText if available.

% -------- manual NV parsing (robust across MATLAB versions)
opts = struct('Mask',[], 'Draw',{{'bbox'}}, 'LineWidth',2, 'Alpha',0.25, ...
              'ShowLabels',true, 'LabelField','confidence', ...
              'Colormap','class', 'FixedColor',[0 255 255], 'SavePath','');
if mod(numel(varargin),2)~=0
    error('Name-value arguments must come in pairs.');
end
for k=1:2:numel(varargin)
    name = varargin{k}; val = varargin{k+1};
    if ~(ischar(name) || (isstring(name)&&isscalar(name)))
        error('Parameter name must be char or string.');
    end
    name = char(name);
    if ~isfield(opts,name), error('Unknown parameter: %s', name); end
    opts.(name) = val;
end

Iout = Irgb;
meta = struct('n_flakes', numel(flakes), 'draw_modes', {opts.Draw}, 'used_cmap', opts.Colormap);
if isempty(flakes), if ~isempty(opts.SavePath), imwrite(Iout, opts.SavePath); end, return; end

H = size(Irgb,1); W = size(Irgb,2);
hasCVT = ~isempty(which('insertShape'));   % Computer Vision Toolbox?

% normalize Draw into a cellstr of lower-case tokens
drawModes = lower(string(opts.Draw));
drawModes = cellstr(drawModes(:));

% colors
if strcmpi(opts.Colormap,'class')
    classes = unique([flakes.class_id]);
    Cmap = classColors(max(1,numel(classes))); % uint8 Kx3
else
    Cmap = uint8(reshape(opts.FixedColor,1,3));
end

for i = 1:numel(flakes)
    f = flakes(i);
    col = selectColor(f.class_id, Cmap, opts.Colormap);

    % ---- bbox
    if any(strcmp(drawModes,'bbox'))
        bb = round(f.bbox_px); % [x y w h]
        if hasCVT
            Iout = insertShape(Iout,'Rectangle',bb,'LineWidth',opts.LineWidth,'Color',double(col));
        else
            Iout = drawRectFallback(Iout, bb, col, opts.LineWidth);
        end
    end

    % ---- rotated (approx from centroid + axes + orientation)
    if any(strcmp(drawModes,'rotated'))
        poly = rotatedRectPoly(f.centroid_px, f.major_axis_px, f.minor_axis_px, f.orientation_deg);
        if hasCVT
            Iout = insertShape(Iout,'Polygon',poly(:).','LineWidth',opts.LineWidth,'Color',double(col));
        else
            Iout = drawPolyFallback(Iout, poly, col, opts.LineWidth);
        end
    end

    % ---- contour
    if any(strcmp(drawModes,'contour'))
        if ~isempty(opts.Mask)
            bb = round(f.bbox_px);
            [x1,x2,y1,y2] = bbox2lims(bb,[H W]);
            sub = opts.Mask(y1:y2, x1:x2);
            if any(sub(:))
                B = bwboundaries(sub,'noholes');
                cx = f.centroid_px(1)-x1+1; cy = f.centroid_px(2)-y1+1;
                b = pickClosestBoundary(B, [cx,cy]);
                if ~isempty(b)
                    poly = [b(:,2)+x1-1, b(:,1)+y1-1];
                    if hasCVT
                        Iout = insertShape(Iout,'Polygon',poly(:).','LineWidth',opts.LineWidth,'Color',double(col));
                    else
                        Iout = drawPolyFallback(Iout, poly, col, opts.LineWidth);
                    end
                end
            else
                % fallback to bbox outline
                bb = round(f.bbox_px);
                if hasCVT
                    Iout = insertShape(Iout,'Rectangle',bb,'LineWidth',opts.LineWidth,'Color',double(col));
                else
                    Iout = drawRectFallback(Iout, bb, col, opts.LineWidth);
                end
            end
        else
            % no mask provided; draw bbox as contour
            bb = round(f.bbox_px);
            if hasCVT
                Iout = insertShape(Iout,'Rectangle',bb,'LineWidth',opts.LineWidth,'Color',double(col));
            else
                Iout = drawRectFallback(Iout, bb, col, opts.LineWidth);
            end
        end
    end

    % ---- label
    if opts.ShowLabels
        labelField = char(opts.LabelField);
        val = NaN;
        if isfield(f,labelField)
            tmp = f.(labelField);
            if isnumeric(tmp) && isscalar(tmp), val = double(tmp); end
        end
        txt = ['#' num2str(f.id) '  ' labelField '=' num2str(val,'%.2f')];
        pos = [round(f.centroid_px(1)+6), max(1, round(f.centroid_px(2)-10))]; % [x y]
        if hasCVT
            Iout = insertText(Iout, pos, txt, 'BoxOpacity',0.6, 'FontSize',14, ...
                              'TextColor','black','BoxColor',double(col));
        else
            Iout = drawDotFallback(Iout, round(f.centroid_px), col);
        end
    end
end

if ~isempty(opts.SavePath)
    outDir = fileparts(char(opts.SavePath));
    if ~exist(outDir,'dir'), mkdir(outDir); end
    imwrite(Iout, opts.SavePath);
end
end

% ================= helpers =================

function col = selectColor(class_id, Cmap, mode)
if strcmpi(mode,'class')
    K = size(Cmap,1);
    if K==0, col = uint8([255 0 0]); return; end
    idx = 1 + mod(max(0,double(class_id)-1), K);
    col = Cmap(idx,:);
else
    col = Cmap(1,:);
end
end

function [x1,x2,y1,y2] = bbox2lims(bb, imsz)
x1 = max(1, round(bb(1)));          y1 = max(1, round(bb(2)));
x2 = min(imsz(2), x1 + max(0, round(bb(3))-1));
y2 = min(imsz(1), y1 + max(0, round(bb(4))-1));
end

function b = pickClosestBoundary(B, c) % c = [cx,cy] in subimage coords
if isempty(B), b = []; return; end
dmin = inf; b = [];
for k=1:numel(B)
    bd = B{k};
    cx = mean(bd(:,2)); cy = mean(bd(:,1));
    d = hypot(cx - c(1), cy - c(2));
    if d<dmin, dmin = d; b = bd; end
end
end

function poly = rotatedRectPoly(centroid, majorAxis, minorAxis, theta_deg)
% approximate rotated rectangle around region stats
w = max(2, minorAxis/2); h = max(2, majorAxis/2);
cx = centroid(1); cy = centroid(2);
theta = deg2rad(-theta_deg); % image y down
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = [ -w -h;  w -h;  w  h; -w  h ]';
rot = R*pts;
poly = [rot(1,:).'+cx, rot(2,:).'+cy];
end

function I2 = drawPolyFallback(I, poly, col, lw)
I2 = I;
edge = polygonEdgeMask(poly, size(I,1), size(I,2), lw);
for c=1:3
    ch = I2(:,:,c);
    ch(edge) = uint8( round(0.7*double(ch(edge)) + 0.3*double(col(c)) ) );
    I2(:,:,c) = ch;
end
end

function I2 = drawRectFallback(I, bb, col, lw)
x1 = bb(1); y1 = bb(2); w = bb(3); h = bb(4);
poly = [x1 y1; x1+w-1 y1; x1+w-1 y1+h-1; x1 y1+h-1];
I2 = drawPolyFallback(I, poly, col, lw);
end

function I2 = drawDotFallback(I, pt, col)
I2 = I;
r = 3; se = strel('disk', r);
M = false(size(I,1), size(I,2));
x = min(max(round(pt(1)),1), size(I,2));
y = min(max(round(pt(2)),1), size(I,1));
M(y,x) = true; M = imdilate(M, se);
for c=1:3
    ch = I2(:,:,c);
    ch(M) = uint8( round(0.2*double(ch(M)) + 0.8*double(col(c)) ) );
    I2(:,:,c) = ch;
end
end

function mask = polygonEdgeMask(poly, H, W, lw)
% rasterize polygon outline with thickness lw
polyX = poly(:,1); polyY = poly(:,2);
maskFill = poly2mask(polyX, polyY, H, W);
edge = bwperim(maskFill);
se = strel('disk', max(1, round(lw/2)));
mask = imdilate(edge, se);
end

function C = classColors(K)
base = [  66 135 245;
         245  66 143;
          60 179 113;
         255 165   0;
         128   0 128;
         220  20  60;
         255 215   0;
          70 130 180;
          46 139  87;
         210 105  30];
if K <= size(base,1)
    C = uint8(base(1:K,:));
else
    reps = ceil(K/size(base,1));
    C = uint8(repmat(base, reps, 1));
    C = C(1:K,:);
end
end
