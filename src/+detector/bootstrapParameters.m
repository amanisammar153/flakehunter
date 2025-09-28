function paths = bootstrapParameters(material, thickness, magnification, baseDir, imageDir, varargin)
%DETECTOR.BOOTSTRAPPARAMETERS  Create minimal Parameters/ from your images.
% Writes:
%   Parameters/GMM_Parameters/<material>_<thickness>.json        (with "1","2",... keys)
%   Parameters/Flatfields/<material>_<thickness>_<mag>x.png      (unity 1x1 flatfield)
%   Parameters/Scan_Magnification/<mag>x.json                    (placeholders)
%
% Name-Value: 'K',2 | 'UsedChannels',"BGR" | 'MaxSamples',20000 | 'RandomSeed',42

ip = inputParser; 
ip.addParameter('K',2); 
ip.addParameter('UsedChannels',"BGR");
ip.addParameter('MaxSamples',20000); 
ip.addParameter('RandomSeed',42); 
ip.parse(varargin{:}); 
opt = ip.Results;
rng(opt.RandomSeed);

baseDir  = char(baseDir); 
imageDir = char(imageDir);
assert(isfolder(baseDir) && isfolder(imageDir), 'Bad baseDir or imageDir.');

% --- make Parameters/ dirs
Pdir = fullfile(baseDir,"Parameters");
gmmDir = fullfile(Pdir,"GMM_Parameters");
ffDir  = fullfile(Pdir,"Flatfields");
magDir = fullfile(Pdir,"Scan_Magnification");
if ~exist(gmmDir,'dir'), mkdir(gmmDir); end
if ~exist(ffDir,'dir'),  mkdir(ffDir);  end
if ~exist(magDir,'dir'), mkdir(magDir); end

% --- unity flatfield (until you capture a real one)
ff = uint8(reshape([128 128 128],[1 1 3]));
ffName = sprintf('%s_%s_%dx.png', lower(material), thickness, magnification);
imwrite(ff, fullfile(ffDir, ffName));

% --- magnification placeholder JSON
mag = struct('view_field_x',0,'view_field_y',0,'x_offset',0,'y_offset',0);
fid=fopen(fullfile(magDir,sprintf('%dx.json',magnification)),'w'); 
fwrite(fid,jsonencode(mag),'char'); fclose(fid);

% --- sample RGB pixels from up to 8 images
ex = ["*.png","*.jpg","*.jpeg","*.tif","*.tiff"];
L = []; for i=1:numel(ex), L=[L; dir(fullfile(imageDir,ex(i)))]; end %#ok<AGROW>
assert(~isempty(L),'No images found in %s', imageDir);
pick = randperm(numel(L), min(numel(L),8));
S = []; target = opt.MaxSamples;

for i = 1:numel(pick)
    I = imread(fullfile(L(pick(i)).folder, L(pick(i)).name));
    if size(I,3)~=3, continue; end
    I = double(I);
    n = numel(I(:,:,1)); take = min(round(target/numel(pick)), n);
    idx = randperm(n, take);
    chIdx = channelIndex(opt.UsedChannels);  % e.g., BGR -> [3 2 1] in MATLAB RGB
    X = [reshape(I(:,:,chIdx(1)),[],1), reshape(I(:,:,chIdx(2)),[],1), reshape(I(:,:,chIdx(3)),[],1)];
    S = [S; X(idx,:)]; %#ok<AGROW>
end
assert(size(S,1) >= 100, 'Not enough RGB samples; add more images or lower MaxSamples.');

% --- fit model: GMM if available, else k-means fallback
useGMM = ~isempty(which('fitgmdist'));
K = opt.K;

if useGMM
    g = fitgmdist(S, K, 'RegularizationValue',1e-6,'Options',statset('MaxIter',500));
    MU = g.mu; SIG = g.Sigma;
else
    warning('Statistics Toolbox not found. Falling back to k-means (approx).');
    [idxK, MU] = kmeans(S, K, 'MaxIter',500, 'Replicates',3); %#ok<ASGLU>
    SIG = zeros(3,3,K);
    for k=1:K
        Sk = S(idxK==k,:); if isempty(Sk), Sk = S; end
        v = var(Sk,0,1);
        SIG(:,:,k) = diag(max(v, 1e-6));
    end
end

% --- heuristic per-channel radius = 2*std per component
RAD = zeros(K,3);
for k=1:K, RAD(k,:) = 2*sqrt(diag(SIG(:,:,k))).'; end

% --- build a *temporary* struct with valid MATLAB field names "x1","x2",...
% then patch JSON keys to "1","2",... after jsonencode.
bgr_idx = channelIndex("BGR"); rgb_idx = channelIndex("RGB");
map = arrayfun(@(ix)find(rgb_idx==ix), bgr_idx);  % [3 2 1]
Jtmp = struct();

for k=1:K
    muJ  = MU(k, map);
    SigJ = SIG(map, map, k);
    radJ = RAD(k, map);
    key  = sprintf('x%d',k); % <-- valid MATLAB field name
    Jtmp.(key) = struct( ...
        'contrast', struct('b',muJ(1),'g',muJ(2),'r',muJ(3)), ...
        'covariance_matrix', reshape(SigJ,1,[]), ...
        'color_radius', struct('b',radJ(1),'g',radJ(2),'r',radJ(3)) );
end

% --- encode and rewrite keys "x1"->"1", "x2"->"2", ...
txt = jsonencode(Jtmp);
txt = regexprep(txt, '"x(\d+)":', '"$1":');

gmmPath = fullfile(gmmDir, sprintf('%s_%s.json', lower(material), thickness));
fid=fopen(gmmPath,'w'); fwrite(fid, txt, 'char'); fclose(fid);

paths = struct('GMM',gmmPath,'Flatfield',fullfile(ffDir,ffName), ...
               'Magnification',fullfile(magDir,sprintf('%dx.json',magnification)));
fprintf('Wrote:\n  %s\n  %s\n  %s\n', paths.GMM, paths.Flatfield, paths.Magnification);
end

function idx = channelIndex(order)
order = upper(string(order)); 
m = containers.Map({'R','G','B'},{1,2,3});  % MATLAB image planes are RGB=1,2,3
letters = char(order);
idx = [m(letters(1)), m(letters(2)), m(letters(3))];
end
