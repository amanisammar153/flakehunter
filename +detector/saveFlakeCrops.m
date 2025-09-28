function files = saveFlakeCrops(Irgb, flakes, outDir, varargin)
ip=inputParser; ip.addParameter('MarginPx',8,@isscalar);
ip.addParameter('Prefix','',@(s)ischar(s)||isstring(s)); ip.parse(varargin{:}); o=ip.Results;
if ~exist(outDir,'dir'), mkdir(outDir); end
H=size(Irgb,1); W=size(Irgb,2); files=strings(0,1);
for i=1:numel(flakes)
    bb = round(flakes(i).bbox_px); x1=max(1,bb(1)-o.MarginPx); y1=max(1,bb(2)-o.MarginPx);
    x2=min(W,bb(1)+bb(3)-1+o.MarginPx); y2=min(H,bb(2)+bb(4)-1+o.MarginPx);
    crop = Irgb(y1:y2, x1:x2, :);
    fp = fullfile(outDir, sprintf('%sflake_%04d_cls%d_conf%03d.png', char(o.Prefix), flakes(i).id, flakes(i).class_id, round(100*flakes(i).confidence)));
    imwrite(crop, fp); files(end+1)=string(fp); %#ok<AGROW>
end
end
