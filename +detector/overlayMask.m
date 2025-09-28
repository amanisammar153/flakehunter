function RGB = overlayMask(I, mask, varargin)
% DETECTOR.OVERLAYMASK  Alpha-blend a binary mask (and its outline) onto RGB.
%
% RGB = detector.overlayMask(I, mask, 'Alpha',0.35, 'EdgeAlpha',1.0, 'EdgeWidth',1, ...
%                            'Color',[1 0 0])
%
% Inputs
%   I      : HxWx3 (uint8/uint16/single/double) base image
%   mask   : HxW logical
% Name–Value
%   'Alpha'     : fill opacity in [0..1] (default 0.35)
%   'EdgeAlpha' : outline opacity in [0..1] (default 1.0)
%   'EdgeWidth' : 0,1,2... pixels (default 1). 0 disables outline.
%   'Color'     : 1x3 RGB in [0..1] (default red [1 0 0])
%
% Output
%   RGB   : same class as I
%
% Notes: toolbox-free (no IPT). Outline via convolution edge detection.

ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Alpha',0.35,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
ip.addParameter('EdgeAlpha',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
ip.addParameter('EdgeWidth',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('Color',[1 0 0],@(v)isnumeric(v)&&numel(v)==3);
ip.parse(varargin{:});
opt = ip.Results;

I0 = I;
if ndims(I0)==2, I0 = repmat(I0,1,1,3); end
if size(I0,3)~=3, error('I must be HxWx3.'); end
cls = class(I0);
I  = double(I0);
mx = double(intmax('uint8'));
if isa(I0,'uint16'), mx = double(intmax('uint16')); end
if isfloat(I0)
    % assume 0..1 or 0..255 range; normalize to 0..1 for blending then restore
    mx = max(1, max(I0(:)));
end

Afill = opt.Alpha;
Aedge = opt.EdgeAlpha;
col   = double(opt.Color(:)).';
col   = max(0,min(1,col));

M = logical(mask);

% outline via 3x3 neighborhood count
edge = false(size(M));
if opt.EdgeWidth>0
    k = true(3); cnt = conv2(double(M), double(k), 'same');
    edge = M & (cnt < 9);
    if opt.EdgeWidth>1
        w = round(opt.EdgeWidth);
        k = ones(2*w+1);
        edge = conv2(double(edge), k, 'same')>0;
    end
end

% composite: first fill, then stronger edge on top
RGBf = I./mx;
for c=1:3
    ch = RGBf(:,:,c);
    ch = ch.*(1 - Afill*M) + (Afill*M)*col(c);
    RGBf(:,:,c) = ch;
end
for c=1:3
    ch = RGBf(:,:,c);
    ch = ch.*(1 - Aedge*edge) + (Aedge*edge)*col(c);
    RGBf(:,:,c) = ch;
end

% restore class
if isfloat(I0)
    RGB = cast(RGBf, cls);
else
    RGB = cast(max(0,min(mx, round(RGBf*mx))), cls);
end
end
