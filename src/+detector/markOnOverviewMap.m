function Iov = markOnOverviewMap(Iov, flakesMM, varargin)
%MARKONOVERVIEWMAP  Plot flake positions (mm) on a stage overview image.
%   Iov = detector.markOnOverviewMap(Iov, flakesMM, 'XRangeMM',105,'YRangeMM',103.333, 'XOffsetMM',0,'YOffsetMM',0, ...)

ip = inputParser;
ip.addParameter('XRangeMM',105,@isscalar);   % stage width covered by the overview image (mm)
ip.addParameter('YRangeMM',103.333,@isscalar);
ip.addParameter('XOffsetMM',0,@isscalar);    % offsets already included in flake positions; here we REMOVE them
ip.addParameter('YOffsetMM',0,@isscalar);
ip.addParameter('RadiusPx',6,@isscalar);
ip.addParameter('LineWidth',2,@isscalar);
ip.addParameter('ShowLabels',true,@islogical);
ip.parse(varargin{:});
opt = ip.Results;

[Hov, Wov, ~] = size(Iov);
hasCVT = ~isempty(which('insertShape'));

% pixels per mm (center-based)
sx = Wov / opt.XRangeMM;
sy = Hov / opt.YRangeMM;

for i=1:numel(flakesMM)
    xm = getfieldOr(flakesMM(i),'flake_position_x_mm',NaN);
    ym = getfieldOr(flakesMM(i),'flake_position_y_mm',NaN);
    if any(isnan([xm,ym])), continue; end

    % remove global offsets so (0,0) maps to image center
    xm = xm - opt.XOffsetMM;
    ym = ym - opt.YOffsetMM;

    % center of overview is (Wov/2, Hov/2); stage y-up → image row decreases
    col = round(Wov/2 + xm * sx);
    row = round(Hov/2 - ym * sy);

    % clamp
    col = min(max(col,1), Wov);
    row = min(max(row,1), Hov);

    % color per class if present, else a default
    colRGB = classColor(uint8(getfieldOr(flakesMM(i),'class_id',1)));

    if hasCVT
        Iov = insertShape(Iov,'FilledCircle',[col row opt.RadiusPx], 'Opacity',0.6, 'Color',double(colRGB));
        if opt.ShowLabels
            txt = sprintf('#%d  conf=%.2f', getfieldOr(flakesMM(i),'id',i), getfieldOr(flakesMM(i),'confidence',NaN));
            Iov = insertText(Iov,[col+8 max(row-10,1)], txt, 'FontSize',14, ...
                             'BoxOpacity',0.6, 'TextColor','black', 'BoxColor',double(colRGB));
        end
    else
        Iov = drawDiskFallback(Iov, row, col, opt.RadiusPx, colRGB);
    end
end
end

% ---- helpers
function v = getfieldOr(S,name,def)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
function col = classColor(k)
base = [66 135 245;245 66 143;60 179 113;255 165 0;128 0 128;220 20 60;255 215 0;70 130 180];
col = uint8(base(1+mod(max(0,double(k)-1), size(base,1)),:));
end
function I = drawDiskFallback(I, r, c, rad, col)
[H,W,~] = size(I);
[xx,yy] = meshgrid(1:W,1:H);
M = (xx-c).^2 + (yy-r).^2 <= rad^2;
I(repmat(M,[1 1 3])) = 255;
end
