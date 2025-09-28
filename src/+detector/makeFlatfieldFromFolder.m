function ffPath = makeFlatfieldFromFolder(imageDir, outPath, varargin)
p=inputParser; p.addParameter('MaxImages',15); p.addParameter('BlurSigmaFrac',1/20);
p.addParameter('TargetMean',128); p.parse(varargin{:}); o=p.Results;
L=[dir(fullfile(imageDir,'*.png'));dir(fullfile(imageDir,'*.jpg'));dir(fullfile(imageDir,'*.jpeg'));dir(fullfile(imageDir,'*.tif'));dir(fullfile(imageDir,'*.tiff'))];
assert(~isempty(L),'No images in %s',imageDir);
idx=round(linspace(1,numel(L),min(numel(L),o.MaxImages)));
I0=im2uint8(imread(fullfile(L(idx(1)).folder,L(idx(1)).name))); [H,W,~]=size(I0);
acc=zeros(H,W,3,'double');
for k=idx(:).'
    I=im2uint8(imread(fullfile(L(k).folder,L(k).name)));
    sig=max(2,round(min(H,W)*o.BlurSigmaFrac));
    for c=1:3, acc(:,:,c)=acc(:,:,c)+imgaussfilt(double(I(:,:,c)),sig); end
end
ff=acc/numel(idx); m=squeeze(mean(mean(ff,1),2)); scale=o.TargetMean./max(m,1);
for c=1:3, ff(:,:,c)=ff(:,:,c)*scale(c); end
ff=uint8(min(max(ff,0),255)); imwrite(ff,outPath); ffPath=outPath;
end
