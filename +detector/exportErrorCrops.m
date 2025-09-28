function OUT = exportErrorCrops(I, predMask, gtMask, outDir, varargin)
% DETECTOR.EXPORTERRORCROPS  Save crops for FP/FN regions with overlays & JSON.
%
% OUT = detector.exportErrorCrops(I, predMask, gtMask, outDir, ...
%          'MarginPx',8, 'MinAreaPx',20, 'ColorFP',[1 0 0], 'ColorFN',[0 0.7 1], ...
%          'Prefix','sample', 'SaveJSON',true)
%
% Inputs
%   I         : HxWx3 image (uint8/uint16/float). Grayscale accepted.
%   predMask  : HxW logical predicted mask
%   gtMask    : HxW logical ground-truth mask
%   outDir    : folder to write crops (created if missing)
%
% Name–Value
%   'MarginPx' : padding around bbox (default 8)
%   'MinAreaPx': skip regions smaller than this (default 20)
%   'ColorFP'  : RGB [0..1] for FP overlay (default red)
%   'ColorFN'  : RGB [0..1] for FN overlay (default cyan)
%   'Prefix'   : filename prefix (default 'sample')
%   'SaveJSON' : also write a JSON manifest (default true)
%
% Output
%   OUT : struct with fields:
%         .fpCount, .fnCount, .dir, .list (table of saved files & boxes)

ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('MarginPx',8,@(x)isscalar(x)&&x>=0);
ip.addParameter('MinAreaPx',20,@(x)isscalar(x)&&x>=0);
ip.addParameter('ColorFP',[1 0 0],@(v)isnumeric(v)&&numel(v)==3);
ip.addParameter('ColorFN',[0 0.7 1],@(v)isnumeric(v)&&numel(v)==3);
ip.addParameter('Prefix','sample',@(s)ischar(s)||isstring(s));
ip.addParameter('SaveJSON',true,@(b)islogical(b)||ismember(b,[0 1]));
ip.parse(varargin{:});
opt = ip.Results;

if ~isfolder(outDir), mkdir(outDir); end
if ndims(I)==2, I = repmat(I,1,1,3); end
I = I(:,:,1:3);

FPmask = logical(predMask & ~gtMask);
FNmask = logical(~predMask & gtMask);

tblRows = [];

% helper to process a mask category
processCat('FP', FPmask, opt.ColorFP);
processCat('FN', FNmask, opt.ColorFN);

T = cell2table(tblRows, 'VariableNames', {'type','file','x0','y0','x1','y1','area_px'});
OUT = struct('fpCount', sum(strcmp(T.type,'FP')), ...
             'fnCount', sum(strcmp(T.type,'FN')), ...
             'dir', outDir, 'list', T);

if opt.SaveJSON
    J = struct('summary', struct('fp',OUT.fpCount,'fn',OUT.fnCount), ...
               'crops', table2struct(T));
    fid = fopen(fullfile(outDir, sprintf('%s_errors.json', char(opt.Prefix))),'w');
    fwrite(fid, jsonencode(J),'char'); fclose(fid);
end

    function processCat(tag, mask, color)
        CC = bwconncomp(mask, 8);
        for k=1:CC.NumObjects
            pix = CC.PixelIdxList{k};
            if numel(pix) < opt.MinAreaPx, continue; end
            [yy,xx] = ind2sub(CC.ImageSize, pix);
            x0 = max(1, min(xx) - opt.MarginPx);
            x1 = min(size(I,2), max(xx) + opt.MarginPx);
            y0 = max(1, min(yy) - opt.MarginPx);
            y1 = min(size(I,1), max(yy) + opt.MarginPx);

            Icrop = I(y0:y1, x0:x1, :);
            Mcrop = false(size(mask));
            Mcrop(pix) = true;
            Mcrop = Mcrop(y0:y1, x0:x1);

            O = detector.overlayMask(Icrop, Mcrop, 'Alpha',0.25, 'EdgeAlpha',1.0, 'EdgeWidth',1, 'Color',color);

            fn = sprintf('%s_%s_%04d_%dx%d.png', char(opt.Prefix), tag, k, x1-x0+1, y1-y0+1);
            fp = fullfile(outDir, fn);
            imwrite(O, fp);
            tblRows(end+1,:) = {tag, fp, x0, y0, x1, y1, numel(pix)}; %#ok<AGROW>
        end
    end
end
